# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later

# Records the state of the git checkout for "achilles --version", "--version-full" and the build report.
#
# Run as a script (cmake -P) on every build by the achilles_git_info target, so the
# information is current after a plain "make", without rerunning CMake. The output is
# only rewritten when something changed, so an unchanged tree does not trigger a rebuild.
#
# Where the information comes from, in order:
#   1. the git checkout itself;
#   2. .achilles-git-info.cmake, a snapshot written into the source tree before packaging
#      (sdists and wheels, see SNAPSHOT below and .github/workflows/pypi-publish.yml);
#   3. CMake/GitArchiveInfo.cmake, filled in by "git archive" (GitHub release tarballs)
#      through the export-subst attribute in .gitattributes.
#
# Inputs:
#   SOURCE_DIR      the Achilles source tree
#   OUTPUT          the C++ file to generate (optional)
#   SNAPSHOT        write the information to this file for packaging (optional)
#   GIT_EXECUTABLE  path to git (may be empty)

cmake_minimum_required(VERSION 3.17)

# safe.directory: inside build containers the checkout often belongs to another user,
# which git refuses to read otherwise
function(run_git out_var)
    execute_process(COMMAND "${GIT_EXECUTABLE}" -c "safe.directory=${SOURCE_DIR}" ${ARGN}
                    WORKING_DIRECTORY "${SOURCE_DIR}"
                    RESULT_VARIABLE result
                    OUTPUT_VARIABLE output
                    ERROR_QUIET
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(NOT result EQUAL 0)
        set(output "")
    endif()
    set(${out_var} "${output}" PARENT_SCOPE)
endfunction()

set(available false)
set(commit "")
set(branch "")
set(describe "")
set(date "")
set(subject "")
set(remote "")
set(dirty false)
set(shallow false)
set(diff_hash "")
set(changed_files "")
set(source "none")

if(GIT_EXECUTABLE)
    run_git(inside rev-parse --is-inside-work-tree)
    # Only describe our own repository, not one that happens to contain the source tree
    run_git(toplevel rev-parse --show-toplevel)
    if(inside STREQUAL "true" AND toplevel)
        get_filename_component(toplevel "${toplevel}" REALPATH)
        get_filename_component(source_dir "${SOURCE_DIR}" REALPATH)
        if(toplevel STREQUAL source_dir)
            set(available true)
        endif()
    endif()
endif()

if(available)
    set(source "git checkout")
    run_git(commit rev-parse HEAD)
    run_git(branch rev-parse --abbrev-ref HEAD)
    # Only release tags (v*); the data release tags (achilles-data-*) are not versions
    run_git(describe describe --tags --match "v*" --always --dirty)
    run_git(date log -1 --format=%cI)
    run_git(subject log -1 --format=%s)
    run_git(is_shallow rev-parse --is-shallow-repository)
    if(is_shallow STREQUAL "true")
        set(shallow true)
    endif()

    # Credentials can be embedded in remote URLs (https://user:token@host/...); drop them
    run_git(remote config --get remote.origin.url)
    string(REGEX REPLACE "://[^/@]+@" "://" remote "${remote}")

    # Tracked files only: untracked files are not part of the build
    run_git(status status --porcelain=v1 --untracked-files=no)
    if(status)
        set(dirty true)
        string(REPLACE "\n" ";" status_lines "${status}")
        foreach(line IN LISTS status_lines)
            # "XY path" or "XY old -> new"
            string(SUBSTRING "${line}" 3 -1 path)
            list(APPEND changed_files "${path}")
        endforeach()
        # Identifies the exact uncommitted changes without including them. Not stripped,
        # so it matches "git diff HEAD --binary | sha256sum" run in the checkout.
        execute_process(COMMAND "${GIT_EXECUTABLE}" -c "safe.directory=${SOURCE_DIR}" diff HEAD --binary
                        WORKING_DIRECTORY "${SOURCE_DIR}"
                        OUTPUT_VARIABLE diff
                        ERROR_QUIET)
        string(SHA256 diff_hash "${diff}")
    endif()
endif()

if(NOT available AND EXISTS "${SOURCE_DIR}/.achilles-git-info.cmake")
    # Snapshot from packaging: sets the same variables as above
    include("${SOURCE_DIR}/.achilles-git-info.cmake")
    set(source "package snapshot")
elseif(NOT available AND EXISTS "${SOURCE_DIR}/CMake/GitArchiveInfo.cmake")
    include("${SOURCE_DIR}/CMake/GitArchiveInfo.cmake")
    # Outside of "git archive" the placeholders are left as they are
    if(NOT ARCHIVE_COMMIT MATCHES "Format:" AND ARCHIVE_COMMIT MATCHES "^[0-9a-f]+$")
        set(available true)
        set(source "release archive")
        set(commit "${ARCHIVE_COMMIT}")
        set(date "${ARCHIVE_DATE}")
        set(subject "${ARCHIVE_SUBJECT}")
        # Older versions of git do not know %(describe) and leave it unexpanded
        if(NOT ARCHIVE_DESCRIBE MATCHES "[%$]")
            set(describe "${ARCHIVE_DESCRIBE}")
        endif()
        # Release tags pointing at the commit, from "HEAD -> main, tag: v1.0.0, origin/main"
        string(REGEX MATCHALL "tag: v[^,]+" tags "${ARCHIVE_REFS}")
        string(REPLACE "tag: " "" tags "${tags}")
        string(REPLACE ";" ", " branch "${tags}")
    endif()
endif()

if(SNAPSHOT)
    # Bracket arguments keep any text in the values as it is
    set(snapshot "# Generated by CMake/WriteGitInfo.cmake for packaging. Do not edit.\n")
    foreach(var available commit branch describe date subject remote dirty shallow diff_hash)
        string(APPEND snapshot "set(${var} [==[${${var}}]==])\n")
    endforeach()
    string(APPEND snapshot "set(changed_files")
    foreach(path IN LISTS changed_files)
        string(APPEND snapshot " [==[${path}]==]")
    endforeach()
    string(APPEND snapshot ")\n")
    file(WRITE "${SNAPSHOT}" "${snapshot}")
    message(STATUS "Wrote git information (${source}) to ${SNAPSHOT}")
endif()

if(NOT OUTPUT)
    return()
endif()

# C++ string literal; the delimiter cannot appear in git output in practice
function(cxx_string out_var value)
    set(${out_var} "R\"git(${value})git\"" PARENT_SCOPE)
endfunction()

cxx_string(commit_s "${commit}")
cxx_string(branch_s "${branch}")
cxx_string(describe_s "${describe}")
cxx_string(date_s "${date}")
cxx_string(subject_s "${subject}")
cxx_string(remote_s "${remote}")
cxx_string(diff_hash_s "${diff_hash}")
cxx_string(source_s "${source}")
set(files_s "")
foreach(path IN LISTS changed_files)
    cxx_string(path_s "${path}")
    string(APPEND files_s "\n        ${path_s},")
endforeach()

set(content "// Generated by CMake/WriteGitInfo.cmake on every build. Do not edit.

#include \"Achilles/GitInfo.hh\"

const achilles::GitInfo &achilles::BuildGitInfo() {
    static const GitInfo info{
        ${available},
        ${commit_s},
        ${branch_s},
        ${describe_s},
        ${date_s},
        ${subject_s},
        ${remote_s},
        ${dirty},
        ${shallow},
        ${diff_hash_s},
        {${files_s}
        },
        ${source_s},
    };
    return info;
}
")

if(EXISTS "${OUTPUT}")
    file(READ "${OUTPUT}" previous)
else()
    set(previous "")
endif()
if(NOT previous STREQUAL content)
    file(WRITE "${OUTPUT}" "${content}")
endif()

# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later

# Filled in by "git archive" through the export-subst attribute (see .gitattributes), which
# is how GitHub creates the source tarballs of a release. In a git checkout the
# placeholders stay as they are and CMake/WriteGitInfo.cmake ignores this file.
set(ARCHIVE_COMMIT [==[$Format:%H$]==])
set(ARCHIVE_DATE [==[$Format:%cI$]==])
set(ARCHIVE_SUBJECT [==[$Format:%s$]==])
set(ARCHIVE_REFS [==[$Format:%D$]==])
set(ARCHIVE_DESCRIBE [==[$Format:%(describe:tags=true,match=v*)$]==])

// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef ACHILLES_TEST_MATCHERS
#define ACHILLES_TEST_MATCHERS

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_templated.hpp"

#include "Achilles/FourVector.hh"
#include "Achilles/ThreeVector.hh"

struct ThreeVectorWithinRelMatcher : Catch::Matchers::MatcherGenericBase {
    ThreeVectorWithinRelMatcher(const achilles::ThreeVector &vector, double err)
        : m_vector{vector}, m_err{err} {}

    bool match(const achilles::ThreeVector &other) const {
        return std::abs(m_vector[0] - other[0]) < m_err &&
               std::abs(m_vector[1] - other[1]) < m_err && std::abs(m_vector[2] - other[2]) < m_err;
    }

    std::string describe() const override { return "WithinRel: " + m_vector.ToString(); }

  private:
    achilles::ThreeVector m_vector;
    double m_err;
};

struct FourVectorWithinRelMatcher : Catch::Matchers::MatcherGenericBase {
    FourVectorWithinRelMatcher(const achilles::FourVector &vector, double err)
        : m_vector{vector}, m_err{err} {}

    bool match(const achilles::FourVector &other) const {
        return std::abs(m_vector[0] - other[0]) < m_err &&
               std::abs(m_vector[1] - other[1]) < m_err &&
               std::abs(m_vector[2] - other[2]) < m_err && std::abs(m_vector[3] - other[3]) < m_err;
    }

    std::string describe() const override { return "WithinRel: " + m_vector.ToString(); }

  private:
    achilles::FourVector m_vector;
    double m_err;
};

inline auto ThreeVectorWithinRel(const achilles::ThreeVector &vec, double err)
    -> ThreeVectorWithinRelMatcher {
    return ThreeVectorWithinRelMatcher(vec, err);
}

inline auto FourVectorWithinRel(const achilles::FourVector &vec, double err)
    -> FourVectorWithinRelMatcher {
    return FourVectorWithinRelMatcher(vec, err);
}

#endif

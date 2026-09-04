// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef ACHILLES_TEST_MATCHERS
#define ACHILLES_TEST_MATCHERS

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_templated.hpp"

#include "Achilles/FourVector.hh"
#include "Achilles/ThreeVector.hh"

template <int D> struct ThreeVectorWithinRelMatcher : Catch::Matchers::MatcherGenericBase {
    ThreeVectorWithinRelMatcher(const achilles::ThreeVectorT<D> &vector, double err)
        : m_vector{vector}, m_err{err} {}

    bool match(const achilles::ThreeVectorT<D> &other) const {
        return std::abs((m_vector[0] - other[0]).native()) < m_err &&
               std::abs((m_vector[1] - other[1]).native()) < m_err &&
               std::abs((m_vector[2] - other[2]).native()) < m_err;
    }

    std::string describe() const override { return "WithinRel: " + m_vector.ToString(); }

  private:
    achilles::ThreeVectorT<D> m_vector;
    double m_err;
};

template <int D> struct FourVectorWithinRelMatcher : Catch::Matchers::MatcherGenericBase {
    FourVectorWithinRelMatcher(const achilles::FourVectorT<D> &vector, double err)
        : m_vector{vector}, m_err{err} {}

    bool match(const achilles::FourVectorT<D> &other) const {
        return std::abs((m_vector[0] - other[0]).native()) < m_err &&
               std::abs((m_vector[1] - other[1]).native()) < m_err &&
               std::abs((m_vector[2] - other[2]).native()) < m_err &&
               std::abs((m_vector[3] - other[3]).native()) < m_err;
    }

    std::string describe() const override { return "WithinRel: " + m_vector.ToString(); }

  private:
    achilles::FourVectorT<D> m_vector;
    double m_err;
};

template <int D>
inline auto ThreeVectorWithinRel(const achilles::ThreeVectorT<D> &vec, double err)
    -> ThreeVectorWithinRelMatcher<D> {
    return ThreeVectorWithinRelMatcher<D>(vec, err);
}

template <int D>
inline auto FourVectorWithinRel(const achilles::FourVectorT<D> &vec, double err)
    -> FourVectorWithinRelMatcher<D> {
    return FourVectorWithinRelMatcher<D>(vec, err);
}

#endif

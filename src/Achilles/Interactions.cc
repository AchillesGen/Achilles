// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#include <map>

#include "Achilles/Event.hh"
#include "Achilles/Interactions.hh"
#include "Achilles/PhysicalUnits.hh"

using namespace achilles;

bool achilles::pid_compare::operator()(const std::pair<PID, PID> &lhs,
                                       const std::pair<PID, PID> &rhs) const {
    PID a1 = lhs.first < lhs.second ? lhs.first : lhs.second;
    PID a2 = lhs.first > lhs.second ? lhs.first : lhs.second;

    PID b1 = rhs.first < rhs.second ? rhs.first : rhs.second;
    PID b2 = rhs.first > rhs.second ? rhs.first : rhs.second;

    return std::tie(a1, a2) < std::tie(b1, b2);
}

// Fit parameters; the values are those of the fit, in the MeV powers noted.
const std::map<std::string, double> HZETRN = {{"a", (5.0_MeV).native()},
                                              {"b", 0.199 / sqrt((1_MeV).native())},
                                              {"c", 0.451 * pow((1_MeV).native(), -0.258)},
                                              {"d", (25.0_MeV).native()},
                                              {"e", (134.0_MeV).native()},
                                              {"f", 1.187 * pow((1_MeV).native(), -0.35)},
                                              {"g", (0.1_MeV).native()},
                                              {"h", (0.282_MeV).native()}};

const std::map<std::string, double> PDG = {{"Zpp", (33.45_mb).in(units::mb)},
                                           {"Zpn", (35.80_mb).in(units::mb)},
                                           {"Y1pp", (42.53_mb).in(units::mb)},
                                           {"Y1pn", (40.15_mb).in(units::mb)},
                                           {"Y2pp", (33.34_mb).in(units::mb)},
                                           {"Y2pn", (30.00_mb).in(units::mb)},
                                           {"B", (0.308_mb).in(units::mb)},
                                           {"s1", 1.0 * pow((1.0_GeV).native(), 2)},
                                           {"s0", pow((5.38_GeV).native(), 2)},
                                           {"n1", 0.458},
                                           {"n2", 0.545}};

const std::map<std::string, double> JWN = {{"gamma", 52.5 * pow((1.0_GeV).native(), 0.16)}, // mb
                                           {"alpha", 0.00369 / (1_MeV).native()},
                                           {"beta", 0.00895741 * pow((1_MeV).native(), -0.8)}};

double Interaction::TotalCrossSection(Event &event, size_t part1, size_t part2) const {
    auto cross_sections = CrossSection(event, part1, part2);
    double total = 0;
    for(const auto &cross_section : cross_sections) total += cross_section.cross_section;
    spdlog::debug("Total cross section for particles {} and {}: {} mb", event.Hadrons()[part1].ID(),
                  event.Hadrons()[part2].ID(), total);
    return total;
}

double Interaction::CrossSectionLab(bool samePID, const double &pLab) const noexcept {
    const double tLab =
        sqrt(pow(pLab, 2) + pow(Constant::mN().native(), 2)) - Constant::mN().native();
    if(samePID) {
        if(pLab < (1.8_GeV).native()) {
            if(tLab >= (25_MeV).native())
                return (1.0 + HZETRN.at("a") / tLab) *
                       (40 + 109.0 * std::cos(HZETRN.at("b") * sqrt(tLab)) *
                                 exp(-HZETRN.at("c") * pow(tLab - HZETRN.at("d"), 0.258)));
            else
                return exp(6.51 * exp(-pow(tLab / HZETRN.at("e"), 0.7)));
        } else if(pLab <= (4.7_GeV).native()) {
            return JWN.at("gamma") / pow(pLab, 0.16);
        } else {
            double ecm2 =
                2 * Constant::mN().native() *
                (Constant::mN().native() + sqrt(pow(pLab, 2) + pow(Constant::mN().native(), 2)));
            return PDG.at("Zpp") + PDG.at("B") * pow(log(ecm2 / PDG.at("s0")), 2) +
                   PDG.at("Y1pp") * pow(PDG.at("s1") / ecm2, PDG.at("n1")) -
                   PDG.at("Y2pp") * pow(PDG.at("s1") / ecm2, PDG.at("n2"));
        }
    } else {
        if(pLab < (0.5_GeV).native()) {
            if(tLab >= (0.1_MeV).native())
                return 38.0 + 12500.0 * exp(-HZETRN.at("f") * pow(tLab - HZETRN.at("g"), 0.35));
            else
                return 26000 * exp(-pow(tLab / HZETRN.at("h"), 0.3));
        } else if(pLab <= (2.0_GeV).native()) {
            return 40 + 10 * cos(JWN.at("alpha") * pLab - 0.943) *
                            exp(-JWN.at("beta") * pow(pLab, 0.8) + 2);
        } else {
            double ecm2 =
                2 * Constant::mN().native() *
                (Constant::mN().native() + sqrt(pow(pLab, 2) + pow(Constant::mN().native(), 2)));
            return PDG.at("Zpn") + PDG.at("B") * pow(log(ecm2 / PDG.at("s0")), 2) +
                   PDG.at("Y1pn") * pow(PDG.at("s1") / ecm2, PDG.at("n1")) -
                   PDG.at("Y2pn") * pow(PDG.at("s1") / ecm2, PDG.at("n2"));
        }
    }
}

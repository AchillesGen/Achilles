#include "nuchic/NucleonPDF.hh"
#include "nuchic/DeuteronDIS.hh"
#include "nuchic/ParticleInfo.hh"
#include "nuchic/Utilities.hh"
#include "nuchic/Vegas.hh"

#include <fstream>

using namespace nuchic;

int main() {
    nuchic::NucleonPDF pdf("MMHT2014lo68cl");
    nuchic::DeuteronDIS deut("data/nk_deut.dat", "MMHT2014lo68cl");
    auto q2vec = nuchic::Linspace(0.5, 2.5, 50);
    constexpr double x = 0.08;
    std::ofstream output_nucleon("F2p_x_008_LO_2.dat");
    for(const auto &q2 : q2vec) {
        auto out = fmt::format("{:8.3e}, {:8.3e}, {:8.3e}\n", q2,
                               pdf.F2(nuchic::PID::proton(), 0.08, q2),
                               pdf.F2(nuchic::PID::proton(), 0.08, q2, 1));
        output_nucleon << out;
    }

    output_nucleon.close();

    return 0;

    YAML::Node config = YAML::LoadFile("dis.yml");
    std::ofstream output("F2d_x_0175.dat");
    output << "q2,val,stat_err,pdf_err+,pdf_err-\n";
    for(const auto &q2 : q2vec) {
        double errup{}, errdw{};
        nuchic::AdaptiveMap map(2);
        nuchic::Vegas integrator(map, config["Initialize"]);
        auto func = [&](const std::vector<double> &ran, const double &wgt) {
            static constexpr double dp = 2;
            static constexpr double dcos = 2;
            double p = dp*ran[0];
            double cos = dcos*ran[1]-1;
            auto results = deut.F2All(x, q2, p, cos);
            for(auto &result : results) result *= dp*dcos; 
            const auto err = pdf.Uncertainty(results); 
            errup += err.errplus * wgt;
            errdw += err.errminus * wgt;
            return results[0];
        };
        integrator(func);
        errup = 0;
        errdw = 0;
        integrator.Set(config["EventGen"]);
        auto result = integrator(func);
        output << fmt::format("{:8.3e}, {:8.3e}, {:8.3e}, {:8.3e}, {:8.3e}\n",
                              q2, result[0]/2, result[1]/2, errdw/2, errup/2);
    }

    output.close();
    return 0;
}

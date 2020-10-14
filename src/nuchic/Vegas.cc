/***************************************************************************
 * Class designed to preform the VEGAS integration algorithm.
 * Contains: Function for VEGAS integration and other necessary functions
 * Inputs: Function to integrate, number of iterations, number of event calls
 * Output: Integrated result, uncertainty
 *
 * The code has been modified to use c++11 threads in order to speed up the
 * code for functions that take a lot of time to evaluate.
 *
 * Inspiration for the VEGAS code comes from the python version written by:
 * Peter Lepage, and can be found at: 
 * https://github.com/gplepage/vegas/blob/master/src/vegas/_vegas.pyx
 **************************************************************************/

#include <algorithm>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <stdexcept>

#include "nuchic/AdaptiveMap.hh"
#include "nuchic/Random.hh"
#include "nuchic/Vegas.hh"

#ifdef USING_MPI
// #include "nuchic/MPI.hh"
#endif

using nuchic::VegasResult;
using nuchic::Vegas;

void VegasResult::Update(double mean_, double var_) {
    sMean.push_back(mean_);
    sVar.push_back(var_);
    if(weighted) {
        meanDivVar += mean_ / var_;
        mean2DivVar += mean_ * mean_ / var_;
        OneDivVar += 1.0 / var_;
        mean = meanDivVar/OneDivVar;
        var = 1.0 / OneDivVar;
    } else {
        sum += mean_;
        sum2 += mean_ * mean_;
        sumVar += var_ + lim::min();
        mean = sum/static_cast<double>(++n);
        var = sumVar/static_cast<double>(n*n);
    }
}

bool VegasResult::Converged(const double& rtol, const double& atol) const {
    return var < pow(atol + rtol * std::abs(mean),2); 
}

double VegasResult::chi2() const {
    if(sMean.size() == 1) return 1;
    if(weighted) return mean2DivVar - meanDivVar * meanDivVar / OneDivVar;
    else return (sum2 - mean*mean * static_cast<double>(n)) * static_cast<double>(n) / var;
}

std::string VegasResult::Summary() const {
    double chi2Dof{};
    if(sMean.size() != 1) {
        chi2Dof = chi2()/dof();
    } else {
        chi2Dof = 0.0;
    }

    return fmt::format("{:3d}   {:^5.3e} +/- {:^5.3e}    {:^5.3e} +/- {:^5.3e}   {:^8.2f}",
            sMean.size(),sMean.back(),std::sqrt(sVar.back()),
            mean,std::sqrt(var), chi2Dof);
}

Vegas::Vegas(nuchic::AdaptiveMap map_, 
             const YAML::Node &args) : map(std::move(map_)) {
    nstrats = 0;
    sumSigF = lim::max();
    lastNEval = 0;
   
    SetDefaults();
    Set(args);
    auto seed = static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    if(args["seed"])
        seed = args["seed"].as<unsigned int>();
#if USING_MPI
    // seed += nuchic::mpi -> Rank();
#endif
    rand = randutils::mt19937_rng(seed);
}

void Vegas::SetDefaults() {
    nitn = nitn_default;
    neval = neval_default;
    maxinc = maxinc_default;
    nCubeBatch = nCubeBatch_default;
    maxCube = maxCube_default;
    maxEvalCube = maxEvalCube_default;
    alpha = alpha_default;
    beta = beta_default;
    adapt = adapt_default;
    adaptErrors = adaptErrors_default;
    rtol = rtol_default;
    atol = atol_default;
    sync = sync_default;
}

void Vegas::Set(const YAML::Node &args) {
    if(args["iterations"])
        nitn = args["iterations"].as<size_t>();
    if(args["evaluations"])
        neval = args["evaluations"].as<size_t>();
    if(args["maxinc"])
        maxinc = args["maxinc"].as<size_t>();
    if(args["nCubeBatch"])
        nCubeBatch = args["nCubeBatch"].as<size_t>();
    if(args["maxCube"])
        maxCube = static_cast<size_t>(args["maxCube"].as<double>());
    if(args["maxEvalCube"])
        maxEvalCube = static_cast<size_t>(args["maxEvalCube"].as<double>());
    if(args["alpha"])
        alpha = args["alpha"].as<double>();
    if(args["beta"])
        beta = args["beta"].as<double>();
    if(args["adapt"])
        adapt = args["adapt"].as<bool>();
    if(args["adaptErrors"])
        adaptErrors = args["adaptErrors"].as<bool>();
    if(args["rtol"])
        rtol = args["rtol"].as<double>();
    if(args["atol"])
        atol = args["atol"].as<double>();
    if(args["sync"])
        sync = args["sync"].as<bool>();

    Setup();
}

void Vegas::Setup() {
    ndims = map.GetDims();
    auto nevalD = static_cast<double>(neval);
    auto ndimsD = static_cast<double>(ndims);
    double nEvalEff = beta > 0 ? (nevalD / 2) : nevalD;
    auto ns = static_cast<size_t>(pow(nEvalEff/2, 1.0/ndimsD)); // Stratifications / axis
    auto ni = static_cast<size_t>(nevalD / 10);                // Increments / axis

    if(ns < 1) ns = 1;
    else if(beta > 0 && static_cast<size_t>(pow(static_cast<double>(ns),ndimsD)) > maxCube) 
        ns = static_cast<size_t>(pow(static_cast<double>(maxCube),1.0/ndimsD));

    if(ni < 1) ni = 1;
    else if(ni > maxinc) ni = maxinc;

    // Want integer number of increments in each stratification and vice versa
    if(ns > ni) {
        if(ns < maxinc) ni = ns;
        else ns = ns / ni * ni;
    } else ni = ni / ns * ns;

    if(adaptErrors) {
        // ni > ns makes no sense with this mode
        if(ni > ns) ni = ns;
    }


    // Rebuild map with correct number of increments
    map.Adapt(0,ni);

    // Determine minimum number of evaluations per cube
    if(nstrats != ns) 
        // Need to recalculate stratification distribution for beta > 0
        sumSigF = lim::max();

    nstrats = ns;
    nCube = static_cast<size_t>(pow(static_cast<double>(nstrats), ndimsD));
    minEvalCube = static_cast<size_t>(nEvalEff / static_cast<double>(nCube));
    if(minEvalCube < 2) minEvalCube = 2;

    if(beta > 0 && sigF.size() != nCube)
        sigF = Batch(nCube,sumSigF/static_cast<double>(nCube));

    nEvalCube = std::vector<size_t>(nCubeBatch,minEvalCube);
}
        
void Vegas::SynchronizeMPI() {
//#ifdef USING_MPI
//    if(nuchic::mpi -> Rank() == 0) 
//         
//#endif
}

nuchic::Batch2D Vegas::GenerateRandom(const std::size_t &dim1, const std::size_t &dim2) {
    nuchic::Batch2D batchResult(dim1);
    for(std::size_t i = 0; i < dim1; ++i) {
        std::vector<double> tmp(dim2);
        rand.generate<std::uniform_real_distribution>(tmp);
        batchResult[i] = tmp;
    }
    return batchResult;
}

void Vegas::RandomBatch(size_t cubeBase, size_t nCubeBatch_, 
        Batch2D& x_, Batch2D& y_, Batch& jac_, BatchInt& cube_) {

    const double dY = 1.0/static_cast<double>(nCube);
    const double nEvalSigF = beta > 0 && sumSigF > 0
        ? static_cast<double>(neval)/2.0/sumSigF : 0.0;
    BatchInt y0(ndims);

    // Determine number of evaluations per cube
    size_t nEvalBatch = 0;
    if(beta > 0) {
        Batch sigF_ = Batch(sigF.begin() + static_cast<int>(cubeBase), sigF.end());
        for(size_t iCube = 0; iCube < nCubeBatch_; ++iCube) {
            nEvalCube[iCube] = size_t(sigF_[iCube]*nEvalSigF);
            if(nEvalCube[iCube] < minEvalCube) 
                nEvalCube[iCube] = minEvalCube;
            else if(nEvalCube[iCube] > maxEvalCube) 
                nEvalCube[iCube] = maxEvalCube;
            nEvalBatch += nEvalCube[iCube];
        }
    } else {
        for(size_t iCube = 0; iCube < nCubeBatch_; ++iCube) nEvalCube[iCube] = minEvalCube;
        nEvalBatch = nCubeBatch_ * minEvalCube;
    }
    lastNEval = nEvalBatch;

    if(y_.size() < nEvalBatch) {
        y_.resize(nEvalBatch);
        x_.resize(nEvalBatch);
        jac_.resize(nEvalBatch);
        fdv2.resize(nEvalBatch);
        if(y_[0].size() != ndims) {
            for(auto &i : y_) i.resize(ndims);
            for(auto &i : x_) i.resize(ndims);
        }
    }
    cube_.resize(nEvalBatch);

    // Generate random points
    Batch2D yran = GenerateRandom(nEvalBatch, ndims);
    size_t iStart = 0;
    for(size_t iCube = 0; iCube < nCubeBatch_; ++iCube) {
        size_t cube = cubeBase + iCube;
        size_t tmpCube = cube;
        for(size_t dim = 0; dim < ndims; ++dim) {
            y0[dim] = tmpCube % nstrats;
            tmpCube = (tmpCube - y0[dim]) / nstrats;
        }
        for(size_t dim = 0; dim < ndims; ++dim) {
            for(size_t i = iStart; i < iStart + nEvalCube[iCube]; ++i) {
                y_[i][dim] = (static_cast<double>(y0[dim]) + yran[i][dim]) 
                           / static_cast<double>(nstrats);
            }
        }
        iStart += nEvalCube[iCube];
    }
    map.Map(y_,x_,jac_);

    // Compute weights and return answers
    iStart = 0;
    for(size_t iCube = 0; iCube < nCubeBatch_; ++iCube) {
        for(size_t i = iStart; i < iStart + nEvalCube[iCube]; ++i) {
            jac_[i] *= dY / static_cast<double>(nEvalCube[iCube]);
            cube_[i] = cubeBase + iCube;
        }
        iStart += nEvalCube[iCube];
    }
}

nuchic::Batch Vegas::MakeBatch(const Func& fcn, Batch2D x_, Batch wgt) {
    Batch res;
    for(size_t i = 0; i < x_.size(); ++i) res.push_back(fcn(x_[i],wgt[i]));
    return res;
}

nuchic::VegasPt Vegas::operator() (const Func& fcn) {
    static const std::string header = fmt::format("{:^3s}   {:^23s}    {:^23s}   {:^8s}\n",
                                                  "itn","integral",
                                                  "wgt average","chi2/dof");
#if USING_MPI
    if(ResBos::mpi -> Rank() == 0) {
#endif
        fmt::print(header);
        std::cout << std::string(header.size(), '-') << std::endl;
#if USING_MPI
    }
#endif

    Batch sigF_ = sigF;

    for(size_t itn = 0; itn < nitn; ++itn) {
        double mean = 0, var = 0, sumSigF_ = 0;

        size_t nCubeBatch_ = std::min(nCubeBatch, nCube);

        for(size_t cubeBase = 0; cubeBase < nCube; cubeBase += nCubeBatch_) {
            if(cubeBase + nCubeBatch_ > nCube) nCubeBatch_ = nCube - cubeBase;

            Batch2D x, y; 
            Batch wgt;
            BatchInt cube;
            RandomBatch(cubeBase,nCubeBatch_,x,y,wgt,cube);

            Batch fdv2_ = fdv2;
            Batch fx = MakeBatch(fcn,x,wgt);

            size_t j = 0;
            for(auto iCube : cube) {
                double sumWF = 0, sumWF2 = 0;
                //int neval = 0;
                while(j < cube.size() && cube[j] == iCube) {
                    double wf = wgt[j]*fx[j];
                    sumWF += wf;
                    sumWF2 += wf*wf;
                    fdv2_[j] = pow(wf*static_cast<double>(nEvalCube[iCube - cube[0]]),2);
                    j++; neval++;
                }
                mean += sumWF;
                auto nevalD = static_cast<double>(neval);
                var += (sumWF2*nevalD - sumWF*sumWF) / (nevalD - 1);
                if(var <= 0) var = mean*mean * lim::min()*10 + lim::min();
                double  sigF2 = std::abs(sumWF2*nevalD - sumWF*sumWF);

                if(beta > 0 && adapt) sumSigF_ += pow(sigF2, beta/2);
                if(adaptErrors && adapt) {
                    fdv2_[j-1] = sigF2;
                    map.AddTrainingData(y[j-1], fdv2_[j-1]);
                }
            }
            if(!adaptErrors && adapt && alpha > 0)
                map.AddTrainingData(y,fdv2_);
        }

        // Accumulate result from this iteration
#if USING_MPI
        double data[2];
        if(adapt) {
            data[0] = mean/var;
            data[1] = 1.0/var;
        } else {
            data[0] = mean;
            data[1] = var;
        }
        ResBos::mpi -> ReduceAll(data, 2, MPI_DOUBLE, MPI_SUM); 
        if(adapt) {
            mean = data[0]/data[1];
            var = 1.0/data[1];
        } else {
            mean = data[0] / ResBos::mpi -> Size();
            var = data[1] / ResBos::mpi -> Size();
        }
        map.MPISync();
#endif
        result.Update(mean,var);

        if(beta > 0 && adapt) sumSigF = sumSigF_;
        if(alpha > 0 && adapt) map.Adapt(alpha);

#if USING_MPI
        if(ResBos::mpi -> Rank() == 0) {
#endif
            std::cout << result.Summary() << std::endl;
#if USING_MPI
        }
#endif

        if(result.Converged(rtol,atol)) break;
    }

#if USING_MPI
    if(ResBos::mpi -> Rank() == 0) {
#endif
        std::cout << std::string(header.size(), '-') << std::endl << std::endl;
#if USING_MPI
    }
#endif

    return result.GetResult();
}

//Below are functions for preforming the modified Kahan Summation, to ensure that the number of threads does
//not effect the result for a fixed random seed
void nuchic::KBNSummation::AddTerm(double value) noexcept { 
    ///Function to add a value to the sum that is being calculated and keep the correction term
    double t = sum + value;
    if(std::abs(sum) >= std::abs(value)){
        correction += ((sum - t) + value);
    } else {
        correction += ((value - t) + sum);
    }
    sum = t;
}



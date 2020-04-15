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
    double chi2Dof;
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
             const std::map<std::string, std::string> &args) : map(std::move(map_)) {
    nstrats = 0;
    sumSigF = lim::max();
    lastNEval = 0;
    
    Set(args);
    auto seed = std::stoul(args.at("seed"));
#if USING_MPI
    // seed += nuchic::mpi -> Rank();
#endif
    rand = new randutils::mt19937_rng(seed);
}

void Vegas::Set(const std::map<std::string,std::string>& args) {
    if(args.find("evaluations") != args.end()) 
        neval = std::stoul(args.at("evaluations"));    
    if(args.find("maxinc") != args.end())
        maxinc = std::stoul(args.at("maxinc"));
    if(args.find("nCubeBatch") != args.end())
        nCubeBatch = std::stoul(args.at("nCubeBatch"));
    if(args.find("maxCube") != args.end())
        maxCube = std::stoul(args.at("maxCube"));
    if(args.find("maxEvalCube") != args.end())
        maxEvalCube = std::stoul(args.at("maxEvalCube"));
    if(args.find("iterations") != args.end())
        nitn = std::stoul(args.at("iterations"));
    if(args.find("alpha") != args.end())
        alpha = std::stod(args.at("alpha"));
    if(args.find("beta") != args.end())
        beta = std::stod(args.at("beta"));
    if(args.find("adapt") != args.end())
        adapt = args.at("adapt") == "true" ? true : false;
    if(args.find("adaptErrors") != args.end())
        adaptErrors = args.at("adaptErrors") == "true" ? true : false;
    if(args.find("vegas_rtol") != args.end())
        rtol = std::stod(args.at("vegas_rtol"));
    if(args.find("vegas_atol") != args.end())
        atol = std::stod(args.at("vegas_atol"));
    if(args.find("sync") != args.end())
        sync = args.at("sync") == "true" ? true : false;

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
        rand -> generate<std::uniform_real_distribution>(tmp);
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

/*
   Vegas::Vegas(const size_t& dim, const int& iseed, const int& itmax, const int& calls, const int& nThreads_) noexcept :
   nDim{dim}, itmx{itmax}, ncall{calls}, nThreads{nThreads_} { ///Vegas initialization function
   for(int iThread = 0; iThread < nThreads; ++iThread) rand.push_back(Random(iseed+iThread));
   init = 0;
   mds = 1;
   nprn = 0;
   alph = 1.5;
   }

   void Vegas::SetRange(std::vector<double> low, std::vector<double> high){ ///Set the range of the variables to integrate over
   if(low.size() != nDim || high.size() != nDim) throw std::runtime_error("Vegas: Integration dimension is inconsistent with ranges given");
   for(size_t i = 0; i < nDim; i++) {
   xl[i] = low[i];
   xu[i] = high[i];
   }
   }

   double Vegas::VEGAS(std::function<double(double*,double)> inFunc) noexcept{ ///Main function for vegas integration
   func = inFunc;
   int k;

   mds = 0;
   if(init <= 0) ///Entry point for a fresh integration
   {
   ndo = 1;
   for(size_t iDim = 0; iDim<nDim; iDim++)
   {
   xi[iDim][0]=1.0;
   }

   if(alph < 0.0001) alph=1.5;
   }

   if(init <= 1) ///Entry point for using the previous grid for a detailed run
   {
   it=0;
   si=0;
   si2=si;
   swgt=si;
   schi=si;
   }

   if(init <= 2){ ///Entry point for adding more points to a previous call to this function
   nd=static_cast<int>(ndmx);
   ng=1;
   if(mds != 0){ ///Preform stratified sampling
   ng=int(pow(ncall/2.0+0.25,1.0/nDim));
   mds=1;
   if(2*ng-static_cast<int>(ndmx) >= 0){
   mds=-1;
   npg=ng/static_cast<int>(ndmx)+1;
   nd=ng/npg;
   ng=npg*nd;
   }
   }
   nCubes = pow(ng,nDim);
   npg=std::max((int)(ncall/nCubes),2);
   calls=(double)npg*(double)nCubes;
   dxg=1.0/ng;
   dv2g=pow(calls*pow(dxg,nDim),2)/npg/npg/(npg-1.0);
   xnd=nd;
   dxg*=xnd;
   xjac=1.0/calls;
   for(size_t iDim = 0; iDim < nDim; iDim++){
   dx[iDim]=xu[iDim]-xl[iDim];
   xjac*=dx[iDim];
   }
   if(nd != ndo){
   for(int i = 0; i < std::max(nd,ndo); i++) r[i]=1.0;
   for(size_t iDim = 0; iDim < nDim; iDim++) rebin(double(ndo)/double(xnd),iDim);
ndo = nd;
}
if(nprn >= 0){
    std::cout << std::left << std::setw(12) << "Iteration" << std::setw(24) <<  "Integral" << std::setw(26) <<  "Accumulated Integral" << std::setw(13) <<
        "Chi**2/it" << std::setw(22) <<  "Time for Iteration" << std::endl;
}
}

//parallelization code
std::vector<std::future<void>> futures(nThreads);
auto boundIntegrate = [this](int calls, int iThread){ return threadIntegrate(calls,iThread); };

for(it=0; it<itmx; it++){
    ti=0.0;
    tsi=ti;
    for(size_t iDim = 0; iDim < nDim; iDim++){
        kg[iDim]=1;
        for(int i = 0; i < static_cast<int>(ndmx); i++){
            D[iDim][i]=0.0;
            Di[iDim][i]=0.0;
        }
    }
    std::chrono::time_point<std::chrono::system_clock> start, end; //used for calculating the time for each iteration
    std::chrono::duration<double> elapsedTime;
    for(;;){
        //reset all the values for each iteration
        fb = 0;
        f2b = 0;
        sfun.Reset();
        sfun2.Reset();
        for(size_t i = 0; i < nDim; i++){
            for(int j = 0; j < static_cast<int>(ndmx); j++){
                d[i][j].Reset(); 
                d[i][j].Reset();
            }
        }

        int threadCalls = npg/nThreads;

        //calculate npg function evaluations over nThreads, and calculate the time it took to preform all the function calls
        start = std::chrono::system_clock::now();
        for(int iThread = 0; iThread < nThreads; iThread++){
            futures[iThread] = std::async(std::launch::async,boundIntegrate,threadCalls,iThread);
        }

        for(int iThread = 0; iThread < nThreads; iThread++){
            futures[iThread].get();
        }
        end = std::chrono::system_clock::now();
        elapsedTime = end-start;

        for(size_t iDim = 0; iDim < nDim; iDim++){
            for(int i = 0; i < static_cast<int>(ndmx); i++){
                D[iDim][i] = d[iDim][i].GetSum();
                Di[iDim][i] = di[iDim][i].GetSum();
            }
        }

        fb = sfun.GetSum();
        f2b = sqrt(sfun2.GetSum()*npg);
        f2b = (f2b-fb)*(f2b+fb);
        if(f2b <= 0.0) f2b= TINY;
        ti += fb;
        tsi += f2b;

        if(mds < 0) {
            for(size_t iDim = 0; iDim < nDim; iDim++) D[iDim][ia[iDim]] += f2b;
        }

        for( k = nDim-1; k >= 0; k--) {
            kg[k] = kg[k]%ng;
            if(++kg[k] != 1) break;
        }
        if(k < 0) break;
    }

    //calculate the value of the integral for this iteration and the accumulated total value
    tsi *= dv2g;
    wgt = 1.0/tsi;
    si += wgt*ti;
    schi += wgt*ti*ti;
    swgt += wgt;
    avgi = si/swgt;
    chi2a = std::max((schi-si*avgi)/(it+0.01),0.0);
    sd = sqrt(1/swgt);
    tsi = sqrt(tsi);
    if(nprn>=0){
        std::cout << std::left << std::setw(12) << it+1 << std::setw(8) << ti << " +/- " << std::setw(11) << tsi << std::setw(8) << avgi << std::setw(5) << " +/-" 
            << std::setw(13) << sd << std::setw(13) << chi2a << std::setw(23) << elapsedTime.count() << std::endl;
        if(nprn>0){
            for(size_t j = 0; j < nDim; j++){
                for(int i = 0; i < nd; i++){
                    std::cout << j <<"\t" << i << "\t" <<  xi[j][i] << "\t" << Di[i][j] << "\t" << D[i][j] << std::endl;
                }
            }
        }
    }

    //Refine the grid
    for(size_t j = 0; j < nDim; j++){
        xo=D[j][0];
        xn=D[j][1];
        D[j][0]=(xo+xn)/2.0;
        dt[j]=D[j][0];
        for(int i = 2; i < nd; i++){
            rc=xo+xn;
            xo=xn;
            xn=D[j][i];
            D[j][i-1]=(rc+xn)/3.0;
            dt[j]+=D[j][i-1];
        }
        D[j][nd-1]=(xn+xo)/2.0;
        dt[j]+=D[j][nd-1];
    }
    for(size_t j = 0; j < nDim; j++){
        rc=0;
        for(int i = 0; i < nd; i++){
            if(D[j][i] < TINY) D[j][i] = TINY;
            r[i]=pow((1.0-D[j][i]/dt[j])/(log(dt[j])-log(D[j][i])),alph);
            rc+=r[i];
        }
        rebin(rc/xnd,j);
    }
}

return avgi;
}

void Vegas::rebin(double ra, int j) noexcept{ ///Refine the grid to improve the accuracy of the integration
    int i;
    int k=0;
    double dr = 0.0, xn = 0.0, xo = 0.0;
    double xin[static_cast<int>(ndmx)];

    for(i=0; i < nd-1; i++){
        while(ra > dr) {
            dr += r[(++k)-1];
        }
        if(k>1) xo=xi[j][k-2];
        xn=xi[j][k-1];
        dr -= ra;
        xin[i]=xn-(xn-xo)*dr/r[k-1];
    }
    for(i=0; i<nd-1;i++) xi[j][i]=xin[i];
    xi[j][nd-1]=1.0;
}

void Vegas::threadIntegrate(int nCubes, int iThread) noexcept{ ///Preforms the part of the function desired to be calculated in parallel
    double lwgt, lxn, lxo, lrc, lx[nDim], lf, lf2;
    int lia[nDim];

    //Loop over the number of calls that need to be made in each thread
    std::vector<double> xArr;
    for(int call = 0; call < nCubes; call++){         
        lwgt = xjac;

        //Obtain random number for this iteration, the random number generator is not thread safe,
        //and therefore needs to be surrounded by a mutex
        xArr = rand[iThread].Get(nDim);

        //Loop over the number of dimensions inside the integral
        for(size_t iDim = 0; iDim < nDim; iDim++){
            lxn = (kg[iDim]-xArr[iDim])*dxg+1.0;
            lia[iDim] = std::max(std::min((int)lxn,static_cast<int>(ndmx)),1);
            if(lia[iDim] > 1) {
                lxo = xi[iDim][lia[iDim]-1] - xi[iDim][lia[iDim]-2];
                lrc = xi[iDim][lia[iDim]-2] + (lxn-lia[iDim])*lxo;
            } else {
                lxo = xi[iDim][lia[iDim]-1];
                lrc = (lxn-lia[iDim])*lxo;
            }
            lx[iDim] = xl[iDim] + lrc*dx[iDim]; 
            lwgt *= lxo*xnd;
        }
        //Calculate the function value and multiple by its weight for the total integral
        lf=lwgt*func(lx,lwgt);
        lf2=lf*lf;

        //Add the results for this function call to the variables containing the total results
        //Note: This has to be done in a thread safe way to ensure the correct results
        std::lock_guard<std::mutex> GridLock(GridRefine);
        sfun.AddTerm(lf);
        sfun2.AddTerm(lf2);
        for(size_t iDim = 0; iDim < nDim; iDim++){
            di[iDim][lia[iDim]].AddTerm(lf);
            if(mds >= 0) {
                d[iDim][lia[iDim]].AddTerm(lf2);
            }
        }
    }
}*/

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



/*
 * Inspiration for the VEGAS code comes from the python version written by:
 * Peter Lepage, and can be found at: 
 * https://github.com/gplepage/vegas/blob/master/src/vegas/_vegas.pyx
 */

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

#include "nuchic/AdaptiveMap.hh"
#include "nuchic/Random.hh"
#include "nuchic/Utilities.hh"

using nuchic::AdaptiveMap;
using nuchic::AdaptiveMap2;

AdaptiveMap::AdaptiveMap(const size_t& ndims_) {
    Vector2D grid_;
    for(size_t dim = 0; dim < ndims_; ++dim) {
        grid_.push_back(std::vector<double>{0,1}); 
    }

    if(grid_[0].size() < 2) 
        throw std::runtime_error("AdaptiveMap: Grid must have at least 2 entries in each dimension.");
    
    InitGrid(grid_,true);

    sumF = Vector2D(grid.size(),std::vector<double>(grid[0].size(),0));
    nF = Vector2D(grid.size(),std::vector<double>(grid[0].size(),0));

//    if(bins[0].size() == 1) MakeUniform(nbins);
//    else Adapt(0,nbins);
}

AdaptiveMap::AdaptiveMap(const Vector2D& grid_, const size_t& nbins_) {
    for(auto dim : grid_) {
        std::sort(dim.begin(), dim.end());
    }

    if(grid_[0].size() < 2) 
        throw std::runtime_error("AdaptiveMap: Grid must have at least 2 entries in each dimension.");
    
    InitGrid(grid_,true);

    sumF = Vector2D(grid.size(),std::vector<double>(grid[0].size(),0));
    nF = Vector2D(grid.size(),std::vector<double>(grid[0].size(),0));

    if(nbins_ != 0 && nbins_ != bins[0].size()) {
        if(bins[0].size() == 1) MakeUniform(nbins_);
        else Adapt(0,nbins_);
    }
}

std::array<double,2> AdaptiveMap::Region(const size_t& dim) const {
    return {grid[dim].front(), grid[dim].back()};
}

std::vector<std::array<double,2>> AdaptiveMap::AllRegion() const {
    std::vector<std::array<double,2>> result;
    for(size_t i = 0; i < grid.size(); ++i) {
        result.push_back(Region(i));
    }
    return result;
}

std::string AdaptiveMap::Settings(const size_t& ngrid) const {
    std::string ans = "";
    if(ngrid > 0) { 
        size_t nskip = nbins/ngrid;
        if(nskip < 1) nskip = 1;
        size_t start = nskip/2;
        for(size_t d = 0; d <  grid.size(); ++d) {
            ans += "grid[" + std::to_string(d) + "] = ";
            for(size_t i = start; i < grid[d].size(); i += nskip) {
                ans += std::to_string(grid[d][i]) + ", ";
            }
        }
    }
    ans.erase(ans.end() - 2);
    return ans;
}

std::vector<double> AdaptiveMap::Random() {
    return std::vector<double>(0); 
}

std::vector<double> AdaptiveMap::operator()(const std::vector<double>& y) {
    jacobian = 1.0;
    std::vector<double> x(y.size(),0);
    for(size_t d = 0; d < grid.size(); ++d) {
        double yNbins = y[d]*static_cast<double>(nbins);
        auto iY = static_cast<size_t>(yNbins);
        double dY = yNbins - static_cast<double>(iY);
        if(iY < nbins) {
            x[d] = grid[d][iY] + bins[d][iY]*dY;
            jacobian *= bins[d][iY]*static_cast<double>(nbins);
        } else {
            x[d] = grid[d][nbins];
            jacobian *= bins[d][nbins-1]*static_cast<double>(nbins);
        }
    }

    return x;
}

void AdaptiveMap::InitGrid(const Vector2D& grid_, const bool& init) {
    ndims = grid_.size();
    nbins = grid_[0].size() - 1;
    grid = grid_;

    if(init || bins.size() != ndims || bins[0].size() != nbins)
        bins = Vector2D(ndims,std::vector<double>(nbins,0));

    for(size_t dim = 0; dim < ndims; ++dim) {
        for(size_t bin = 0; bin < nbins; ++bin) {
            bins[dim][bin] = grid[dim][bin+1] - grid[dim][bin];
        }
    }
}

void AdaptiveMap::MakeUniform(const size_t& nbins_) {
    if(nbins_ == 0) nbins = bins[0].size();
    else nbins = nbins_;
    if(nbins < 1) 
        throw std::runtime_error("AdaptiveMap: Requesting " + std::to_string(nbins) 
                + ", but minimum is 1");

    Vector2D newGrid(ndims,std::vector<double>(nbins+1,0));
    for(size_t dim = 0; dim < ndims; ++dim) {
        newGrid[dim] = Linspace(grid[dim].front(), grid[dim].back(), nbins_+1);
    }

    InitGrid(newGrid);
}

void AdaptiveMap::AddTrainingData(const Vector2D& y, const Vector& f) {
    if(sumF[0].size() != nbins) {
        for(auto &i : sumF) i.resize(nbins);
        for(auto &i : nF) i.resize(nbins);
    }
    for(size_t dim = 0; dim < ndims; ++dim) {
        for(size_t i = 0; i < y.size(); ++i) {
            auto iY = static_cast<size_t>(y[i][dim]*static_cast<double>(nbins));
            sumF[dim][iY] += std::abs(f[i]);
            nF[dim][iY] += 1;
        }
    }
}

void AdaptiveMap::AddTrainingData(const Vector& y, const double& f) {
    if(sumF[0].size() != nbins) {
        for(auto &i : sumF) i.resize(nbins);
        for(auto &i : nF) i.resize(nbins);
    }
    for(size_t dim = 0; dim < ndims; ++dim) {
        auto iY = static_cast<size_t>(y[dim]*static_cast<double>(nbins));
        sumF[dim][iY] += std::abs(f);
        nF[dim][iY] += 1;
    }
}

void AdaptiveMap::Map(Vector2D& y, Vector2D& x, Vector& jac) {
    for(size_t i = 0; i < x.size(); ++i) {
        jac[i] = 1.0;
        for(size_t dim = 0; dim < ndims; ++dim) {
            double yInc = y[i][dim] * static_cast<double>(nbins);
            auto iY = static_cast<size_t>(std::floor(yInc));
            double dy = yInc - static_cast<double>(iY);
            if(iY < nbins) {
                x[i][dim] = grid[dim][iY] + bins[dim][iY]*dy;
                jac[i] *= bins[dim][iY] * static_cast<double>(nbins);
            } else {
                x[i][dim] = grid[dim][nbins];
                jac[i] *= bins[dim][nbins-1]*static_cast<double>(nbins);
            }
        }
    }
}

void AdaptiveMap::MPISync() {
#if USING_MPI
    size_t size = ndims*nbins;
    double *vals = new double[2*size];
    for(size_t dim = 0; dim < ndims; ++dim) {
        for(size_t bin = 0; bin < nbins; ++bin) {
            vals[dim*nbins+bin] = sumF[dim][bin];
            vals[dim*nbins+bin+size] = nF[dim][bin];
        }
    }
    ResBos::mpi -> ReduceAll(vals, 2*size, MPI_DOUBLE, MPI_SUM);
    for(size_t dim = 0; dim < ndims; ++dim) {
        for(size_t bin = 0; bin < nbins; ++bin) {
            sumF[dim][bin] = vals[dim*nbins+bin];
            nF[dim][bin] = int(vals[dim*nbins+bin+size]);
        }
    }
    delete vals;
#endif
}

void AdaptiveMap::Adapt(const double& alpha, const size_t& nbins_) {
    size_t oldBins = nbins;
    size_t newBins;
    if(nbins_ == 0) newBins = oldBins;
    else newBins = nbins_;

    if(newBins < 1) 
        throw std::runtime_error("AdaptiveMap: Requesting " + std::to_string(newBins) 
                + ", but minimum is 1");
   
    Vector2D newGrid(ndims,std::vector<double>(newBins+1,0));
    if(newBins == 1) {
        for(size_t dim = 0; dim < ndims; ++dim) {
            newGrid[dim][0] = grid[dim].front();
            newGrid[dim][1] = grid[dim].back();
        }
        InitGrid(newGrid);
        return;
    }

    std::vector<double> averageF(oldBins,1), tmpF;
    if(alpha > 0 && oldBins > 1) {
        tmpF.resize(oldBins,0);
    }
    for(size_t dim = 0; dim < ndims; ++dim) {
        if(!sumF.empty() && alpha != 0) 
            for(size_t i = 0; i < oldBins; ++i) {
                if(nF[dim][i] != 0) averageF[i] = sumF[dim][i]/nF[dim][i];
                else averageF[i] = 0;
            }

        if(alpha > 0 && oldBins > 1) {
            tmpF[0] = (3*averageF[0]+averageF[1])/4;
            tmpF[tmpF.size()-1] = (3*averageF[averageF.size()-1]+averageF[averageF.size()-2])/4;
            double totSumF = tmpF[0] + tmpF[tmpF.size()-1];
            for(size_t i = 1; i < oldBins-1; ++i) {
                tmpF[i] = (6*averageF[i] + averageF[i-1] + averageF[i+1]) / 8;
                totSumF += tmpF[i];
            }

            if(totSumF > 0) 
                for(size_t i = 0; i < oldBins; ++i) averageF[i] = tmpF[i] / totSumF + lim::min();
            else
                for(size_t i = 0; i < oldBins; ++i) averageF[i] = lim::min();

            for(size_t i = 0; i < oldBins; ++i) 
                averageF[i] = pow(-(1-averageF[i])/log(averageF[i]),alpha);
        }

        newGrid[dim][0] = grid[dim][0];
        newGrid[dim][newGrid[0].size()-1] = grid[dim][grid[0].size()-1];
        auto j = static_cast<size_t>(-1);
        double accumF = 0, fBins = 0;
        for(size_t i = 0; i < oldBins; ++i) fBins += averageF[i];
        fBins /= static_cast<double>(newBins);
        for(size_t i = 1; i < newBins; ++i) {
            while(accumF < fBins) {
                j++;
                if(j < oldBins) accumF += averageF[j];
                else break;
            }
            if(accumF >= fBins) {
                accumF -= fBins;
                newGrid[dim][i] = grid[dim][j+1]-(accumF/averageF[j])*bins[dim][j];
                continue;
            }
            break;
        }
    }

    InitGrid(newGrid);
    sumF = Vector2D(grid.size(),std::vector<double>(grid[0].size(),0));
    nF = Vector2D(grid.size(),std::vector<double>(grid[0].size(),0));
}

bool AdaptiveMap2::Deserialize(std::istream &in) {
    in >> m_bins >> m_dims;
    m_hist.resize((m_bins + 1) * m_dims);

    for(auto &x : m_hist) in >> x;

    return true;
}

bool AdaptiveMap2::Serialize(std::ostream &out) const {
    out << m_bins << ' ' << m_dims;

    for(const auto &x : m_hist) {
        out << ' ' << std::scientific
            << std::setprecision(std::numeric_limits<double>::max_digits10 - 1) << x;
    }

    return true;
}

size_t AdaptiveMap2::FindBin(size_t dim, double x) const {
    const auto edges = Edges(dim);
    auto it = std::lower_bound(edges.begin(), edges.end(), x);
    return static_cast<size_t>(std::distance(edges.begin(), it))-1;
}

double AdaptiveMap2::operator()(std::vector<double> &rans) {
    double jacobian = 1.0;
    for(std::size_t i = 0; i < m_dims; ++i) {
        const auto position = rans[i] * static_cast<double>(m_bins);
        const auto index = static_cast<size_t>(position);
        const auto loc = position - static_cast<double>(index);
        const double left = lower_edge(i, index);
        const double size = width(i, index);

        // Calculate inverse CDF
        rans[i] = left + loc * size;

        jacobian *= size * static_cast<double>(m_bins);
    }

    return jacobian;
}

double AdaptiveMap2::GenerateWeight(const std::vector<double> &rans) const {
    double jacobian = 1.0;
    for(std::size_t i = 0; i < m_dims; ++i) {
        const auto index = FindBin(i, rans[i]);
        jacobian *= width(i, index) * static_cast<double>(m_bins);
    }
    return jacobian;
}

void AdaptiveMap2::Adapt(const double &alpha, const std::vector<double> &data) {
    std::vector<double> tmp(m_bins);
    std::vector<double> new_hist(m_hist.size());

    for(size_t i = 0; i < m_dims; ++i) {
        // Load data into tmp
        tmp.assign(data.begin() + static_cast<int>(i * m_bins),
                   data.begin() + static_cast<int>((i + 1) * m_bins));

        // Smooth data by averaging over neighbors
        double previous = tmp[0];
        double current = tmp[1];
        tmp[0] = (previous + current) / 2;
        double norm = tmp[0];

        for(size_t bin = 1; bin < m_bins - 1; ++bin) {
            const double sum = previous + current;
            previous = current;
            current = tmp[bin + 1];
            tmp[bin] = (sum + current) / 3;
            norm += tmp[bin];
        }
        tmp.back() = (previous + current) / 2;
        norm += tmp.back();

        // If norm is zero, then there is no data to adjust
        if(norm == 0) continue;

        // Compute the importance factor
        double average_bin = 0;
        for(size_t bin = 0; bin < m_bins; ++bin) {
            if(tmp[bin] != 0) {
                const double r = tmp[bin] / norm;
                const double fac = pow((r - 1.0) / log(r), alpha);
                average_bin += fac;
                tmp[bin] = fac;
            }
        }
        average_bin /= static_cast<double>(m_bins);

        double cbin = 0;
        size_t ibin = 0;

        // Adjust boundaries
        for(size_t nbin = 1; nbin < m_bins; ++nbin) {
            for(; cbin < average_bin; ++ibin) cbin += tmp[ibin]; 

            const double prev = lower_edge(i, ibin-1);
            const double curr = lower_edge(i, ibin);

            cbin -= average_bin;
            const double delta = (curr - prev) * cbin;
            new_hist[i * (m_bins + 1) + nbin] = curr - delta / tmp[ibin - 1];
        }
        new_hist[i * (m_bins + 1) + m_bins] = 1.0;
    }

    m_hist = new_hist;
}

void AdaptiveMap2::Split(nuchic::AdaptiveMapSplit split) {
    size_t nsplit{};
    if(split == AdaptiveMapSplit::half) {
        nsplit = 2;
    } else if(split == AdaptiveMapSplit::third) {
        nsplit = 3;
    } else if(split == AdaptiveMapSplit::quarter) {
        nsplit = 4; 
    } else {
        throw std::logic_error("Invalid histogram splitting");
    }

    // Create temporary histogram to store new histogram
    std::vector<double> hist(m_dims*(nsplit*m_bins + 1));

    // Split the histogram into new binnings
    for(size_t d = 0; d < m_dims; ++d) {
        for(size_t bin = 0; bin < m_bins; ++bin) {
            const double curr = lower_edge(d, bin);
            const double size = width(d, bin);

            for(size_t i = 0; i < nsplit; ++i) {
                const size_t idx = d*(nsplit*m_bins+1) + bin*nsplit + i;
                hist[idx] = curr + static_cast<double>(i)*size/static_cast<double>(nsplit);
            }
        }
        // Add the endpoint
        hist[(d+1)*(nsplit*m_bins)+d] = 1.0;
    }

    // Store the new histogram information
    m_bins = nsplit*m_bins;
    m_hist = hist;
}

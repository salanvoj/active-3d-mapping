#include <boost/chrono.hpp>
#include <cmath>
#include <mex.h>
#include <omp.h>
#include <queue>
#include <vector>

namespace
{

class Timer
{
private:
    typedef boost::chrono::high_resolution_clock Clock;
    typedef Clock::time_point Time;
    typedef boost::chrono::duration<double> Duration;
    Time start;
public:
    Timer(): start(Clock::now()) {}
    Timer(const Time &s): start(s) {}
    Timer(const Timer &s): start(s.start) {}
    void reset() { start = Clock::now(); }
    double secondsElapsed() const
    {
        return boost::chrono::duration_cast<Duration>(Clock::now() - start).count();
    }
};

typedef double D;
typedef size_t I;
typedef std::vector<D> VecD;
typedef std::vector<I> VecI;
typedef std::vector<VecD> VecVecD;
typedef std::vector<VecI> VecVecI;

void fillRange(const I n, VecI &vec)
{
    vec.resize(n);
    for (I i = 0; i < n; ++i)
        vec[i] = i;
}

template<typename C>
class IndexValueComp
{
public:
    typedef typename C::value_type T;
    IndexValueComp(const C &cont): cont(cont) {}
    bool operator()(const T &a, const T &b) {
        assert(a < cont.size());
        assert(b < cont.size());
        if (std::isnan(cont[a]) && !std::isnan(cont[b]))
            return true;
        return cont[a] < cont[b];
    }
private:
    const C &cont;
};
typedef IndexValueComp<VecD> IndexValueCompVecD;

template<typename It>
class IterValueComp
{
public:
    typedef typename std::iterator_traits<It>::difference_type I;
    IterValueComp(It iter): iter(iter) {}
    bool operator()(const I &a, const I &b) {
        if (std::isnan(*(iter + a)) && !std::isnan(*(iter + b)))
            return true;
        return *(iter + a) < *(iter + b);
    }
private:
    const It iter;
};
typedef IterValueComp<VecD::iterator> IterValueCompD;

typedef std::priority_queue<I, VecI, IterValueCompD> PriorityQueueI;


template<typename It0, typename It1, typename It2>
double computeViewImprovement(It0 v_it, It0 v_it_end, It1 p_it, It2 b)
{
    typename std::iterator_traits<It1>::value_type d = 0.0;
    for (; v_it < v_it_end; ++v_it, ++p_it)
        d += b[*v_it] * (*p_it);
    return d;
}

template<typename It0, typename It1, typename It2>
double updateVoxelCosts(It0 v_it, It0 v_it_end, It1 p_it, It2 b)
{
    for (; v_it < v_it_end; ++v_it, ++p_it)
        b[*v_it] *= (1.0 - *p_it);
}

///
/// @brief planGreedyPrioritized
///
/// @param v        Voxel indices corresponding to p
/// @param p        Probability p[j][i] of voxel v[j][i] being visible in view j.
/// @param b_in     Initial voxel costs
/// @param M        Number of voxels
/// @param L        Number of positions
/// @param K        Number of views to select per position
/// @param J        Indices of selected views (index past the last ray denotes not selected)
/// @param b        Final voxel costs
///
void planGreedyPrioritized(const VecVecI &v,
                           const VecVecD &p,
                           const D *b_in, const I M,
                           const I L, const I K,
                           D *J, D *b_sum, D *b)
{
    // Check and preprocess arguments.
    assert(v.size() == p.size());
    const I NL = v.size();          // Number of rays in total
    const I N = NL / L;             // Number of rays per position
    assert(K <= N);
    std::fill_n(J, L * K, NL);      // Start with no selected rays.
    std::copy(b_in, b_in + M, b);
    // Initialize local variables.
    VecI q(NL);                     // Prioritized ray permutation
    fillRange(NL, q);               // Initialized from 0 to (NL-1).
    VecD d(NL);                     // Ray improvements
    Timer t;
    #pragma omp parallel for schedule(runtime)
    for (size_t j = 0; j < v.size(); ++j)
    {
        d[j] = computeViewImprovement(v[j].begin(), v[j].end(), p[j].begin(), b);
    }
//    std::cout << "plan_rays_greedy: Initialization of ray gains: " << t.secondsElapsed() << "s." << std::endl;
    t.reset();
    size_t n_deltas = v.size();     // Number of recomputed improvements
    PriorityQueueI queue(IterValueCompD(d.begin()), q);
    VecI J_sum(L, 0);               // Number of rays selected at each per position
    size_t k = 0;                   // Number of selected rays
    while (k < L * K) {
        I j = queue.top();
//        std::cout << "Testing ray " << j << " with value " << d[j] << " as " << k << "th ray." << std::endl;
        queue.pop();
        if (d[j] == 0)              // No improvement possible
            break;
        I l = j / N;
        if (J_sum[l] == K)          // Discard rays at saturated positions.
            continue;
        n_deltas++;
        // Recompute ray improvement (from a previous iteration).
        d[j] = computeViewImprovement(v[j].begin(), v[j].end(), p[j].begin(), b);
        if (d[j] < d[queue.top()])  // Others may be better.
        {
            queue.push(j);
            continue;
        }
        // As good as it gets, ray selected.
        k++;
        *J++ = j;
        J_sum[l]++;
        // Update voxel costs.
        updateVoxelCosts(v[j].begin(), v[j].end(), p[j].begin(), b);
    }
//    std::cout << "plan_rays_greedy: Ray selection and gain recomputation: " << t.secondsElapsed() << "s." << std::endl;
//    std::cout << "plan_rays_greedy: " << k << " selected rays, " << n_deltas << " ray gains (re)computed." << std::endl;
}


} // namespace

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    if (nrhs < 4)
        mexErrMsgTxt("Invalid number of arguments.");
    if (!mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]))
        mexErrMsgTxt("Invalid argument for A.");
    if (!(mxGetM(prhs[1]) == 1 || mxGetN(prhs[1]) == 2) || !mxIsDouble(prhs[1]))
        mexErrMsgTxt("Invalid argument for c.");
    if (mxGetM(prhs[0]) != mxGetNumberOfElements(prhs[1]))
        mexErrMsgTxt("Incompatible sizes of A and c.");
    if (!mxIsScalar(prhs[2]) || !mxIsDouble(prhs[2]))
        mexErrMsgTxt("Invalid argument for n_pos.");
    if (!mxIsScalar(prhs[3]) || !mxIsDouble(prhs[3]))
        mexErrMsgTxt("Invalid argument for n_sel.");

    const size_t n_pos_ray = mxGetN(prhs[0]);
    const double *pr = mxGetPr(prhs[0]);
    const size_t *ir = mxGetIr(prhs[0]);
    const size_t *jc = mxGetJc(prhs[0]);
    VecVecI ray_v(n_pos_ray);
    VecVecD ray_p(n_pos_ray);
    for (size_t i = 0; i < n_pos_ray; ++i)
    {
        const size_t n_ray_vox = jc[i + 1] - jc[i];
        ray_v[i].resize(n_ray_vox);
        ray_p[i].resize(n_ray_vox);
        for (size_t j = 0; j < n_ray_vox; ++j)
        {
            ray_v[i][j] = ir[jc[i] + j];
            ray_p[i][j] = pr[jc[i] + j];
        }
    }

    const size_t n_vox = mxGetNumberOfElements(prhs[1]);
    const double *v_val = mxGetPr(prhs[1]);

    const size_t n_pos = mxGetScalar(prhs[2]);
    const size_t n_sel = mxGetScalar(prhs[3]);

    // Output (n_pos * n_sel) indices instead of prob. mask.
    // plhs[0] = mxCreateDoubleMatrix(n_pos_ray, 1, mxREAL);
    plhs[0] = mxCreateDoubleMatrix(n_pos * n_sel, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(n_vox, 1, mxREAL);
//     double *sel_prob = mxGetPr(plhs[0]);
    double *sel_ind = mxGetPr(plhs[0]);
    double *val = mxGetPr(plhs[1]);
    double *vox_p = mxGetPr(plhs[2]);

    planGreedyPrioritized(ray_v, ray_p, v_val, n_vox, n_pos, n_sel, sel_ind, val, vox_p);
}

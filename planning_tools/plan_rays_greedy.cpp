
#include <algorithm>
#include <boost/chrono.hpp>
#include <cassert>
#include <limits>
#include <cmath>
#include <matrix.h>
#include <mex.h>
#include <numeric>
#include <omp.h>
#include <vector>

#include <iostream>
#include <fstream>

//#define SAVE_TIMINGS

namespace {

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
        if (!std::isnan(cont[a]) && std::isnan(cont[b]))
            return true;
        return (cont[a] > cont[b]);
    }
private:
    const C &cont;
};

//void planGreedyPriority(const VecVecI &ray_v,
//                        const VecVecD &ray_p,
//                        const double *v_val, const size_t n_vox,
//                        const size_t n_pos, const size_t n_sel,
//                        double *sel_prob, double *val, double *vox_p)
//{
//    // Timeings and delta updates.
//    Timer t;
//    double t_sorting = 0.0;
//    size_t n_selected = 0;
//    size_t n_deltas = 0;
//    std::ofstream t_sorting_out("/tmp/active_mapping/t_sorting.txt", std::ios::out | std::ios::app);
//    std::ofstream n_deltas_out("/tmp/active_mapping/n_deltas.txt", std::ios::out | std::ios::app);
//    std::ofstream n_selected_out("/tmp/active_mapping/n_selected.txt", std::ios::out | std::ios::app);

//    // Number of positions times number of rays per position.
//    const I n_pos_ray = ray_v.size();
//    const I n_ray = n_pos_ray / n_pos;
//    // Prob. of voxels not being visible.
//    VecD v_neg(n_vox, 1.0);
//    // Number of rays selected at each position.
//    VecI sel(n_pos, 0);
//    // Prob. of rays being selected.
//    std::fill_n(sel_prob, n_pos_ray, 0.0);
//    // Current ray value, start with infinity so that all values are recomputed in the first iteration.
//    VecD r_val(n_pos_ray, std::numeric_limits<D>::infinity());
//    // Value ordering, starting with order imposed by index.
//    VecI r_ord(n_pos_ray);
//    fillRange(n_pos_ray, r_ord);

//    I pos_left = n_pos;
//    while (pos_left > 0)
//    {
//        VecI::iterator it;
//        for (it = r_ord.begin(); it != r_ord.end(); ++it)
//        {
//            n_deltas++;
//            assert(*it < r_val.size());
//            assert(*it < n_vox);
//            // Zero first, before optional skipping.
//            r_val[*it] = 0.0;
//            // Do not count rays already used.
//            if (sel_prob[*it] > 0.0)
//                continue;
//            // Do not count positions already saturated.
//            if (sel[*it / n_ray] == n_sel)
//                continue;
//            VecI::const_iterator v_it;
//            VecD::const_iterator p_it;
//            for (v_it = ray_v[*it].begin(), p_it = ray_p[*it].begin(); v_it != ray_v[*it].end(); ++v_it, ++p_it)
//                r_val[*it] += v_val[*v_it] * (*p_it) * v_neg[*v_it];
//            // Halt if the recomputed value is better than a previous value of the next ray.
//            if (it + 1 != r_ord.end() && r_val[*it] > r_val[*(it + 1)])
//            {
//                // Move to the first sorted element.
//                ++it;
//                break;
//            }
//        }
//        // Sort indices using updated values.
//        // std::sort(r_ord.begin(), r_ord.end(), IndexValueComp<VecD>(r_val));
//        // Use the fact that the tail is already sorted.
//        t.reset();
//        std::sort(r_ord.begin(), it, IndexValueComp<VecD>(r_val));
//        std::inplace_merge(r_ord.begin(), it, r_ord.end(), IndexValueComp<VecD>(r_val));
//        t_sorting += t.secondsElapsed();
//        I i_max = r_ord.front();
//        I pos = i_max / n_ray;
//        // Halt if all voxels have been covered.
//        if (r_val[i_max] == 0.0 || sel[pos] == n_sel)  // Value is zero for saturated positions.
//            break;
//        assert(r_val[i_max] != 0.0);
//        assert(sel[pos] < n_sel);
//        if (++sel[pos] == n_sel)
//            --pos_left;
//        sel_prob[i_max] = 1.0;
//        n_selected++;
//        // Update probabilities of the voxels not being visible.
//        VecI::const_iterator v_it;
//        VecD::const_iterator p_it;
//        for (v_it = ray_v[i_max].begin(), p_it = ray_p[i_max].begin();
//             v_it != ray_v[i_max].end(); ++v_it, ++p_it)
//            v_neg[*v_it] *= (1.0 - *p_it);
//    }
//    // TODO: Complete plan if (pos_left > 0)?
//    *val = 0.0;
//    for (I i = 0; i < n_vox; ++i)
//    {
//        vox_p[i] = 1.0 - v_neg[i];
//        *val += v_val[i] * vox_p[i];
//    }

//    //    mexPrintf("Delta updates: %lu\n", n_deltas);
//    //    mexPrintf("Average delta updates: %.3f\n", double(n_deltas) / n_selected);

//    t_sorting_out << t_sorting << std::endl;
//    n_selected_out << n_selected << std::endl;
//    n_deltas_out << n_deltas << std::endl;

//    t_sorting_out.close();
//    n_selected_out.close();
//    n_deltas_out.close();
//}


///
/// \brief planGreedyPrioritized
///
/// \param v        Voxel indices corresponding to p
/// \param p_vis    Probability p[j][i] of voxel i being visible in ray j
/// \param b_in     Initial voxel costs
/// \param M        Number of voxels
/// \param L        Number of positions
/// \param K        Number of selected rays per position
/// \param J_ind    Selected rays as indicator vector
/// \param b_sum    Sum of voxel costs b
/// \param b        Final voxel costs
///
void planGreedyPrioritized(const VecVecI &v,
                           const VecVecD &p_vis,
                           const double *b_in, const size_t M,
                           const size_t L, const size_t K,
                           double *J_ind, double *b_sum, double *b)
{
#ifdef SAVE_TIMINGS
    // Timings and delta updates.
    Timer t;
    double t_sorting = 0.0;
    size_t n_selected = 0;
    size_t n_deltas = 0;
    std::ofstream t_sorting_out("/tmp/active_mapping/t_sorting.txt", std::ios::out | std::ios::app);
    std::ofstream n_deltas_out("/tmp/active_mapping/n_deltas.txt", std::ios::out | std::ios::app);
    std::ofstream n_selected_out("/tmp/active_mapping/n_selected.txt", std::ios::out | std::ios::app);
#endif
//    std::ofstream log_out("/tmp/active_mapping/plan_rays_greedy.log", std::ios::out | std::ios::app);

    const I NL = v.size();          // Number of rays in total
    const I N = NL / L;             // Number of rays per position
    VecI J_sum(L, 0);               // Number of rays selected at each per position
    std::fill_n(J_ind, NL, 0.0);    // Start with no selected rays.
    VecD d(NL, std::numeric_limits<double>::infinity());  // Force recompute.
    VecI q(NL);                     // Prioritized ray permutation
    fillRange(NL, q);               // Initialized from 0 to (NL-1).
    std::copy(b_in, b_in + M, b);

    size_t k = 0;
    for (; k < L * K; ++k)
    {
        VecI::iterator it;
        for (it = q.begin() + k; it != q.end(); ++it)
        {
#ifdef SAVE_TIMINGS
            n_deltas++;
#endif
            d[*it] = 0.0;
            if (J_sum[*it / N] == K)  // Skip the rays at saturated positions.
                continue;
            VecI::const_iterator v_it;
            VecD::const_iterator p_it;
            for (v_it = v[*it].begin(), p_it = p_vis[*it].begin(); v_it != v[*it].end(); ++v_it, ++p_it)
                d[*it] += b[*v_it] * (*p_it);
            if (it + 1 != q.end() && d[*it] > d[*(it + 1)])
                break;
        }

#ifdef SAVE_TIMINGS
        t.reset();
#endif
        size_t n = std::distance(q.begin(), it);
        std::sort(q.begin() + k, q.begin() + n, IndexValueComp<VecD>(d));
        std::inplace_merge(q.begin() + k, q.begin() + n, q.end(), IndexValueComp<VecD>(d));
#ifdef SAVE_TIMINGS
        t_sorting += t.secondsElapsed();
#endif
        I j_k = q[k];
        if (d[j_k] == 0.0)      // Break if no further improvement can be made.
            break;
        I l_k = j_k / N;
        assert(J_sum[l_k] < K); // Saturated positions were skipped with d = 0.
        J_sum[l_k]++;
        J_ind[j_k] = 1.0;
        // Update current voxel costs.
        VecI::const_iterator v_it;
        VecD::const_iterator p_it;
        for (v_it = v[j_k].begin(), p_it = p_vis[j_k].begin(); v_it != v[j_k].end(); ++v_it, ++p_it)
            b[*v_it] *= (1.0 - *p_it);
    }
    if (k < L * K)
        std::cout << "Only " << k << " rays planned from " << L * K << " allowed." << std::endl;
    *b_sum = std::accumulate(b, b + M, 0.0);
#ifdef SAVE_TIMINGS
    t_sorting_out << t_sorting << std::endl;
    n_selected_out << k << std::endl;
    n_deltas_out << n_deltas << std::endl;
    t_sorting_out.close();
    n_selected_out.close();
    n_deltas_out.close();
#endif
}

///
/// \brief planGreedyPrioritized
///
/// \param v        Voxel indices corresponding to p
/// \param p_vis    Probability p[j][i] of voxel i being visible in ray j
/// \param b_in     Initial voxel costs
/// \param M        Number of voxels
/// \param L        Number of positions
/// \param K        Number of selected rays per position
/// \param J_ind    Selected rays as indicator vector
/// \param b_sum    Sum of voxel costs b
/// \param b        Final voxel costs
///
void planGreedyPrioritizedRemoving(const VecVecI &v,
                                   const VecVecD &p_vis,
                                   const double *b_in, const size_t M,
                                   const size_t L, const size_t K,
                                   double *J_ind, double *b_sum, double *b)
{
    std::cout << "planGreedyPrioritizedRemoving" << std::endl;
#ifdef SAVE_TIMINGS
    // Timings and delta updates.
    Timer t;
    double t_sorting = 0.0;
    size_t n_selected = 0;
    size_t n_deltas = 0;
    std::ofstream t_sorting_out("/tmp/active_mapping/t_sorting.txt", std::ios::out | std::ios::app);
    std::ofstream n_deltas_out("/tmp/active_mapping/n_deltas.txt", std::ios::out | std::ios::app);
    std::ofstream n_selected_out("/tmp/active_mapping/n_selected.txt", std::ios::out | std::ios::app);
#endif
//    std::ofstream log_out("/tmp/active_mapping/plan_rays_greedy.log", std::ios::out | std::ios::app);

    const I NL = v.size();          // Number of rays in total
    const I N = NL / L;             // Number of rays per position
    VecI J_sum(L, 0);               // Number of rays selected at each per position
    std::fill_n(J_ind, NL, 0.0);    // Start with no selected rays.
    VecD d(NL, std::numeric_limits<double>::infinity()); // Current ray delta
    VecI q(NL);                     // Prioritized ray permutation
    fillRange(NL, q);               // Initialized from 0 to (NL-1).
    std::copy(b_in, b_in + M, b);

    while (!q.empty())
    {
        VecI::iterator it;
        for (it = q.begin(); it != q.end(); ++it)
        {
#ifdef SAVE_TIMINGS
            n_deltas++;
#endif
            d[*it] = 0.0;
            VecI::const_iterator v_it;
            VecD::const_iterator p_it;
            for (v_it = v[*it].begin(), p_it = p_vis[*it].begin(); v_it != v[*it].end(); ++v_it, ++p_it)
                d[*it] += b[*v_it] * (*p_it);
            if (it + 1 != q.end() && d[*it] >= d[*(it + 1)])
                break;
        }

#ifdef SAVE_TIMINGS
        t.reset();
#endif
        size_t n = std::distance(q.begin(), it);
        // The boundary ray will be sorted no matter in which subsequence it remains.
        std::sort(q.begin(), q.begin() + n, IndexValueComp<VecD>(d));
        std::inplace_merge(q.begin(), q.begin() + n, q.end(), IndexValueComp<VecD>(d));
#ifdef SAVE_TIMINGS
        t_sorting += t.secondsElapsed();
#endif
        I j_best = q.front();          // Add the best ray
        I l_best = j_best / N;
        J_ind[j_best] = 1.0;
        VecI::const_iterator v_it;  // Update voxel costs
        VecD::const_iterator p_it;
        for (v_it = v[j_best].begin(), p_it = p_vis[j_best].begin(); v_it != v[j_best].end(); ++v_it, ++p_it)
            b[*v_it] *= (1.0 - *p_it);
        J_sum[l_best]++;
        if (J_sum[l_best] == K)        // Close position
        {
            VecI q_new;
            q_new.reserve(q.size());
            for (VecI::const_iterator it = q.begin(); it != q.end(); ++it)
                if (*it / N != l_best)
                    q_new.push_back(*it);
            q.swap(q_new);
        }
        else
            q.erase(q.begin());
    }
    *b_sum = std::accumulate(b, b + M, 0.0);
#ifdef SAVE_TIMINGS
    t_sorting_out << t_sorting << std::endl;
    n_selected_out << K * L << std::endl;
    n_deltas_out << n_deltas << std::endl;
    t_sorting_out.close();
    n_selected_out.close();
    n_deltas_out.close();
#endif
}

void planGreedy(const VecVecI &ray_v,
                const VecVecD &ray_p,
                const double *v_val, const size_t n_vox,
                const size_t n_pos, const size_t n_sel,
                double *sel_prob, double *val, double *vox_p)
{
#ifdef SAVE_TIMINGS
    // Timings and delta updates.
    Timer t;
    size_t n_selected = 0;
    size_t n_deltas = 0;
    std::ofstream n_deltas_out("/tmp/active_mapping/n_deltas.txt", std::ios::out | std::ios::app);
    std::ofstream n_selected_out("/tmp/active_mapping/n_selected.txt", std::ios::out | std::ios::app);
#endif

    // Number of positions times number of rays per position.
    const I n_pos_ray = ray_v.size();
    const I n_ray = n_pos_ray / n_pos;
    // Prob. of voxels not being visible.
    VecD v_neg(n_vox, 1.0);
    // Number of rays selected at each position.
    VecI sel(n_pos, 0);
    // Prob. of rays being selected.
    std::fill_n(sel_prob, n_pos_ray, 0.0);
    VecD r_val(n_pos_ray, 0.0);

    I pos_left = n_pos;
    while (pos_left > 0)
    {
#ifndef SAVE_TIMINGS
#pragma omp parallel for schedule(runtime)
#endif
        for (I i = 0; i < n_pos_ray; ++i)
        {
#ifdef SAVE_TIMINGS
            n_deltas++;
#endif
            // Zero first, before optional skipping.
            r_val[i] = 0.0;
            // Do not count rays already used.
            if (sel_prob[i] > 0.0)
                continue;
            // Do not count positions already saturated.
            if (sel[i / n_ray] == n_sel)
                continue;
            VecI::const_iterator v_it;
            VecD::const_iterator p_it;
            for (v_it = ray_v[i].begin(), p_it = ray_p[i].begin(); v_it != ray_v[i].end(); ++v_it, ++p_it)
                r_val[i] += v_val[*v_it] * (*p_it) * v_neg[*v_it];
        }
        I i_max = std::distance(r_val.begin(), std::max_element(r_val.begin(), r_val.end()));
        I pos = i_max / n_ray;
        // Halt if all voxels have been covered.
        if (r_val[i_max] == 0.0 || sel[pos] == n_sel)
            break;
        assert(r_val[i_max] != 0.0);
        assert(sel[pos] < n_sel);
        if (++sel[pos] == n_sel)
            --pos_left;
        sel_prob[i_max] = 1.0;
#ifdef SAVE_TIMINGS
        n_selected++;
#endif
        // Update probabilities of the voxels not being visible.
        VecI::const_iterator v_it;
        VecD::const_iterator p_it;
        for (v_it = ray_v[i_max].begin(), p_it = ray_p[i_max].begin(); v_it != ray_v[i_max].end(); ++v_it, ++p_it)
            v_neg[*v_it] *= (1.0 - *p_it);
    }
    // TODO: Complete plan if (pos_left > 0)?
    *val = 0.0;
    for (I i = 0; i < n_vox; ++i)
    {
        vox_p[i] = 1.0 - v_neg[i];
        *val += v_val[i] * vox_p[i];
    }
#ifdef SAVE_TIMINGS
    n_selected_out << n_selected << std::endl;
    n_deltas_out << n_deltas << std::endl;
    n_selected_out.close();
    n_deltas_out.close();
#endif
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
//    if (nrhs >= 5 && (!mxIsScalar(prhs[4]) || !mxIsLogical(prhs[4])))
    if (nrhs >= 5 && !mxIsScalar(prhs[4]))
        mexErrMsgTxt("Invalid argument for prioritized.");

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

    // Default to prioritized search.
//    bool priority = (nrhs < 5) || mxGetScalar(prhs[4]);
    size_t priority = (nrhs < 5) ? 2 : mxGetScalar(prhs[4]);

    plhs[0] = mxCreateDoubleMatrix(n_pos_ray, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(n_vox, 1, mxREAL);
    double *sel_prob = mxGetPr(plhs[0]);
    double *val = mxGetPr(plhs[1]);
    double *vox_p = mxGetPr(plhs[2]);

    if (priority >= 2)
        planGreedyPrioritizedRemoving(ray_v, ray_p, v_val, n_vox, n_pos, n_sel, sel_prob, val, vox_p);
    else if (priority == 1)
        planGreedyPrioritized(ray_v, ray_p, v_val, n_vox, n_pos, n_sel, sel_prob, val, vox_p);
    else
        planGreedy(ray_v, ray_p, v_val, n_vox, n_pos, n_sel, sel_prob, val, vox_p);
}

#pragma once

#include <fstream>
#include <Kokkos_Random.hpp>

#include "SimInfo.h"
#include "BoundaryConditions.h"

namespace fv2d {

namespace {

  using RandomPool = Kokkos::Random_XorShift64_Pool<>;

  /**
   * @brief Sod Shock tube aligned along the X axis
   */
  KOKKOS_INLINE_FUNCTION
  void initSodX(Array Q, int i, int j, const Params &params) {
    if (getPos(params, i, j)[IX] <= 0.5) {
      Q(j, i, IR) = 1.0;
      Q(j, i, IP) = 1.0;
      Q(j, i, IU) = 0.0;
    }
    else {
      Q(j, i, IR) = 0.125;
      Q(j, i, IP) = 0.1;
      Q(j, i, IU) = 0.0;
    }
  }

  /**
   * @brief Sod Shock tube aligned along the Y axis
   */
  KOKKOS_INLINE_FUNCTION
  void initSodY(Array Q, int i, int j, const Params &params) {
    if (getPos(params, i, j)[IY] <= 0.5) {
      Q(j, i, IR) = 1.0;
      Q(j, i, IP) = 1.0;
      Q(j, i, IU) = 0.0;
    }
    else {
      Q(j, i, IR) = 0.125;
      Q(j, i, IP) = 0.1;
      Q(j, i, IU) = 0.0;
    }
  }

  /**
   * @brief Sedov blast initial conditions
   */
  KOKKOS_INLINE_FUNCTION
  void initBlast(Array Q, int i, int j, const Params &params) {
    real_t xmid = 0.5 * (params.xmin+params.xmax);
    real_t ymid = 0.5 * (params.ymin+params.ymax);

    Pos pos = getPos(params, i, j);
    real_t x = pos[IX];
    real_t y = pos[IY];

    real_t xr = xmid - x;
    real_t yr = ymid - y;
    real_t r = sqrt(xr*xr+yr*yr);

    if (r < 0.2) {
      Q(j, i, IR) = 1.0;
      Q(j, i, IU) = 0.0;
      Q(j, i, IP) = 10.0;
    }
    else {
      Q(j, i, IR) = 1.2;
      Q(j, i, IU) = 0.0;
      Q(j, i, IP) = 0.1;
    }
  }

  /**
   * @brief Stratified convection based on Hurlburt et al 1984
   */
  KOKKOS_INLINE_FUNCTION
  void initH84(Array Q, int i, int j, const Params &params, const RandomPool &random_pool) {
    Pos pos = getPos(params, i, j);
    real_t x = pos[IX];
    real_t y = pos[IY];

    real_t T = y;
    real_t rho = pow(y, params.m1);
    real_t prs = pow(y, params.m1+1.0); 

    auto generator = random_pool.get_state();
    real_t pert = params.h84_pert * (generator.drand(-0.5, 0.5));
    random_pool.free_state(generator);

    Q(j, i, IR) = rho;
    Q(j, i, IU) = 0.0;
    Q(j, i, IV) = pert;
    Q(j, i, IP) = prs;
  }

  /**
   * @brief Stratified convection based on Cattaneo et al. 1991
   */
  KOKKOS_INLINE_FUNCTION
  void initC91(Array Q, int i, int j, const Params &params, const RandomPool &random_pool) {
    Pos pos = getPos(params, i, j);
    real_t x = pos[IX];
    real_t y = pos[IY];

    real_t T = (1.0 + params.theta1*y);
    real_t rho = pow(T, params.m1);
    real_t prs = pow(T, params.m1+1.0);

    auto generator = random_pool.get_state();
    real_t pert = params.c91_pert * (generator.drand(-0.5, 0.5));
    random_pool.free_state(generator);

    prs = prs * (1.0 + pert);

    Q(j, i, IR) = rho;
    Q(j, i, IU) = 0.0;
    Q(j, i, IV) = 0.0;
    Q(j, i, IP) = prs;
  }

  /**
   * @brief Simple diffusion test with a structure being advected on the grid
   */
  KOKKOS_INLINE_FUNCTION
  void initDiffusion(Array Q, int i, int j, const Params &params) {
    real_t xmid = 0.5 * (params.xmin+params.xmax);
    real_t ymid = 0.5 * (params.ymin+params.ymax);

    Pos pos = getPos(params, i, j);

    real_t x0 = (pos[IX]-xmid);
    real_t y0 = (pos[IY]-ymid);

    real_t r = sqrt(x0*x0+y0*y0);

    if (r < 0.2) 
      Q(j, i, IR) = 1.0;
    else
      Q(j, i, IR) = 0.1;

    Q(j, i, IP) = 1.0;
    Q(j, i, IU) = 1.0;
    Q(j, i, IV) = 1.0;
  }

  /**
   * @brief Rayleigh-Taylor instability setup
   */
  KOKKOS_INLINE_FUNCTION
  void initRayleighTaylor(Array Q, int i, int j, const Params &params) {
    real_t ymid = 0.5*(params.ymin + params.ymax);

    Pos pos = getPos(params, i, j);
    real_t x = pos[IX];
    real_t y = pos[IY];

    const real_t P0 = 2.5;

    if (y < ymid) {
      Q(j, i, IR) = 1.0;
      Q(j, i, IU) = 0.0;
      Q(j, i, IP) = P0 + 0.1 * params.g * y;
    }
    else {
      Q(j, i, IR) = 2.0;
      Q(j, i, IU) = 0.0;
      Q(j, i, IP) = P0 + 0.1 * params.g * y;
    }

    if (y > -1.0/3.0 && y < 1.0/3.0)
      Q(j, i, IV) = 0.01 * (1.0 + cos(4*M_PI*x)) * (1 + cos(3.0*M_PI*y))/4.0;
  }

  /**
   * @brief Reading a spline from the disk
   **/
   void initProfile(Array Q, const Params &params) {
    // Reading input file
    std::string filename = params.init_filename;
    std::vector<real_t> y, rho, u, v, p;
    std::ifstream f_in(filename);

    while (f_in.good()) {
      real_t y_, rho_, u_, v_, p_;
      f_in >> y_ >> rho_ >> u_ >> v_ >> p_;
      if (f_in.good()) {
        y.push_back(y_);
        rho.push_back(rho_);
        u.push_back(u_);
        v.push_back(v_);
        p.push_back(p_);
      }
    }
    f_in.close();
    
    // Copying profile on GPU
    size_t N = y.size();
    Kokkos::View<real_t**> profile("profile", N, 5);
    auto profile_host = Kokkos::create_mirror_view(profile);

    std::cout << "Profile read from " << filename << " has " << N << " points" << std::endl;

    for (size_t i=0; i < N; ++i) {
      profile_host(i, 0) = y[i];
      profile_host(i, 1) = rho[i];
      profile_host(i, 2) = u[i];
      profile_host(i, 3) = v[i];
      profile_host(i, 4) = p[i];
    }

    Kokkos::deep_copy(profile, profile_host);

    // Initializing domain
    Kokkos::parallel_for("Initialization from profile",
                         params.range_dom,
                         KOKKOS_LAMBDA(const int i, const int j) {
                          auto pos = getPos(params, i, j);
                          real_t y = pos[IY];
                    
                          // Finding current cell position in profile.
                          // Could be optimized if dy in the profile is fixed
                          int iy = 0;
                          real_t prof_y = profile(iy, 0);
                          constexpr real_t eps = 1.0e-5;
                          while (prof_y-eps < y && Kokkos::abs(prof_y - y) > eps) {
                            iy++;
                            prof_y = profile(iy, 0);
                          }
                    
                          // Linear interpolation
                          real_t fy = (y - profile(iy-1, 0)) / (profile(iy, 0) - profile(iy-1, 0));
                          for (int ivar=1; ivar < 5; ++ivar)
                            Q(j, i, ivar-1) = profile(iy-1, ivar) * (1.0 - fy) + profile(iy, ivar) * fy;
                        
                          // Todo : Edge case extrapolation
                         });
  }
}



/**
 * @brief Enum describing the type of initialization possible
 */
enum InitType {
  SOD_X,
  SOD_Y,
  BLAST,
  RAYLEIGH_TAYLOR,
  DIFFUSION,
  H84,
  C91,
  HSE,
  PROFILE
};

struct InitFunctor {
private:
  Params params;
  InitType init_type;
public:
  InitFunctor(Params &params)
    : params(params) {
    std::map<std::string, InitType> init_map {
      {"sod_x", SOD_X},
      {"sod_y", SOD_Y},
      {"blast", BLAST},
      {"rayleigh-taylor", RAYLEIGH_TAYLOR},
      {"diffusion", DIFFUSION},
      {"H84", H84},
      {"C91", C91},
      {"profile", PROFILE}
    };

    if (init_map.count(params.problem) == 0)
      throw std::runtime_error("Error unknown problem " + params.problem);

    init_type = init_map[params.problem];
  };
  ~InitFunctor() = default;

  void init(Array &Q) {
    auto init_type = this->init_type;
    auto params = this->params;

    RandomPool random_pool(params.seed);

    // Filling active domain ...
    Kokkos::parallel_for( "Initialization", 
                          params.range_dom, 
                          KOKKOS_LAMBDA(const int i, const int j) {
                            switch(init_type) {
                              case SOD_X:           initSodX(Q, i, j, params); break;
                              case SOD_Y:           initSodY(Q, i, j, params); break;
                              case BLAST:           initBlast(Q, i, j, params); break;
                              case DIFFUSION:       initDiffusion(Q, i, j, params); break;
                              case RAYLEIGH_TAYLOR: initRayleighTaylor(Q, i, j, params); break;
                              case H84:             initH84(Q, i, j, params, random_pool); break;
                              case C91:             initC91(Q, i, j, params, random_pool); break;
                              default: break;
                            }
                          });

    // If filling is via a spline
    if (init_type == PROFILE)
      initProfile(Q, params);
  
    // ... and boundaries
    BoundaryManager bc(params);
    bc.fillBoundaries(Q);
  }
};



}
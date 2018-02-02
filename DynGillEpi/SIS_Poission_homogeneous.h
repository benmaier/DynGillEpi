#ifndef __SIS_POISS_HOMO_H__
#define __SIS_POISS_HOMO_H__
#include <Utilities.h>

SI_result
    SIS_Poisson_homogeneous(size_t N,
                            CONTACTS_LIST contactListList,
                            double infection_rate_per_dt,
                            double recovery_rate_per_dt
                            size_t T_simulation,
                            size_t output_time_resolution,
                            size_t number_of_simulations = 1,
                            size_t initial_number_of_infected = 1,
                            size_t seed = 0,
                            size_t t_infection_start = 0,
                            bool verbose = false
            );

#endif

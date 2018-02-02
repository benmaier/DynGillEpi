/* Simulates independent realizations of a homogeneous SIR process on
a temporal network given on the form: (t i j) with one triple per line
and t,i, and j separated by tabs. ("\t").

Should be compiled with g++ including the boost library and using the
option -O2, i.e., as
g++ SIS-Poisson-homogeneous.cpp -o SIS -O2 -I<path of boost>

With the program compiled as SIS, it is called from the shell as:
./SIS <data> dt beta mu T_simulation ensembleSize outputTimeResolution
where:
<data> - path of text file containing contact data (temporal network);
dt - time-resolution of recorded contact data (time-step length);
beta - probability per time-step of infection when in contact with an
    infected node;
mu - probability per time-step of recovery for an infected node;
T_simulation - length of simulation in number of time-steps;
ensembleSize - number of independent realizations of the SIR process;
outputTimeResolution - time-resolution of the average number of infected
    and recovered nodes that the program gives as output.

The program gives as output two text files containing the average
number of infected nodes as function of time and a histogram of the number of
recovered nodes at t=T_simulation-1.*/
//======================================================================
// Libraries
//======================================================================
#include <Utilities.h>



//======================================================================
// Main:
//========================= =============================================
SI_result
    SIS_Poisson_homogeneous(size_t N,
                            CONTACTS_LIST contactListList,
                            double infection_rate_per_dt,
                            double recovery_rate_per_dt
                            size_t T_simulation,
                            size_t output_time_resolution,
                            size_t number_of_simulations,
                            size_t initial_number_of_infected,
                            size_t seed,
                            size_t t_infection_start,
                            bool verbose
            )
{
    // Set parameter values as specified:
    double beta = infection_rate;
    double mu = recovery_rate;
    COUNTER ensembleSize = number_of_simulations; //ensemble size (number of realizations)
    COUNTER outputTimeResolution = output_time_resolution; //output time-resolution

    //-------------------------------------------------------------------------------------
    // Define variables:
    //-------------------------------------------------------------------------------------
    NODES infected; //list of infected nodes
    double Mu; //cumulative recovery rate
    BOOLS isInfected; //list which nodes are infected
    COUNTER I; //number of infected nodes
    COUNTER SI; //number of susceptible nodes in contact with infectious nodes
    NODES si_s; //list of susceptible nodes in contact with infected nodes
    double Beta; //total infection rate
    double Lambda; //cumulative transition rate
    double xi;
    COUNTER t; //time countet
    COUNTER t_infectionStart = t_infection_start; //starting time of infection
    NODE root; //root node of infection
    double tau; //renormalized waiting time until next event
    NODE i,j; //nodes
    CONTACTS::iterator contact_iterator; //iterator over contacts
    CONTACTS_LIST::iterator contactList_iterator; //iterator over list of contacts
    double r_transitionType; //random variable for choosing which transition happens
    COUNTER m; //transition process
    COUNTER n; //time counter
    NODES::iterator node_iterator; //iterator over list of nodes
    NODES::iterator last; //iterator for use when generating unique list of new infected nodes
    // Containers for output data:
    vector < vector < size_t > > sumI_t(number_of_simulations,T_simulation/outputTimeResolution); //list of number of infected nodes in each recorded frame
    vector < vector < size_t > > sumSI_t(number_of_simulations,T_simulation/outputTimeResolution); //list of number of infected nodes in each recorded frame
    vector < size_t > hist_I(number_of_simulations); //histogram of R values after I=0
    // Random number generators:
    //
    if (seed==0)
    {
        seed = time(nullptr);
    }
    ENG generator(seed);
    DIST_REAL rand(0,1); //random float on [0,1[
    DIST_EXP randexp(1.0); //random exponentially distributed float

    //-------------------------------------------------------------------------------------
    // Simulate:
    //-------------------------------------------------------------------------------------
    std::clock_t start = std::clock();     //timer
    COUNTER stopped=0; //counter of number of simulations that stopped (I=0) during T_simu
    for(int q=0; q<ensembleSize; q++)
    {
        if (verbose)
            std::cout << q << "/" << ensembleSize << std::endl; //print realization # to screen

        // Choose at random infectious root nodes and run SIS process starting from roots:
        vector < size_t > infected_nodes_for_shuffling;
        for(size_t n=0; n<N;++n)
            infected_nodes_for_shuffling.push_back(n);
        choose_random_unique(
                                infected_nodes_for_shuffling.begin(),
                                infected_nodes_for_shuffling.end(),
                                initial_number_of_infected,
                                generator,
                                rand
                            );
        infected.clear();
        isInfected.assign(N,false);
        for(size_t n=0; n<initial_number_of_infected; ++n)
        {
            infected.push_back(infected_nodes_for_shuffling[n]);
            isInfected[n] = true;
        }

        I = initial_number_of_infected;
        Mu = mu*initial_number_of_infected;

       // First waiting time:
        tau = randexp(1);
        // set simulation time to zero:
        t = 0;

        //--- Loop over list of contact lists: ---
        while(I>0 && t<T_simulation) //loop until either I=0 or t>=T_simu
        {
            for(contactList_iterator=contactListList.begin()+t_infectionStart; contactList_iterator!=contactListList.end(); contactList_iterator++)
            {
                // Create list of susceptible nodes in contact with infected nodes:
                si_s.clear();
                for(contact_iterator=(*contactList_iterator).begin(); contact_iterator!=(*contactList_iterator).end(); contact_iterator++)
                {
                    i=(*contact_iterator).first;
                    j=(*contact_iterator).second;
                    if(isInfected[i])
                    {
                        if(!isInfected[j])
                        {
                            si_s.push_back(j);
                        }
                    }
                    else
                    {
                        if(isInfected[j])
                        {
                            si_s.push_back(i);
                        }
                    }
                }
                SI=si_s.size(); //number of possible S->I transitions
                Beta=(double)SI*beta; //cumulative infection rate
                Lambda=Beta+Mu; //cumulative transition rate

                // Check if transition takes place during time-step:
                if(tau>=Lambda) //no transition takes place
                {
                    tau-=Lambda;
                }
                else //at least one transition takes place
                {
                    xi=1.; //fraction of time-step left before transition
                    // Sampling step:
                    while(tau<xi*Lambda) //repeat if next tau is smaller than ~ Lambda-tau
                    {
                        xi-=tau/Lambda; //fraction of time-step left after transition
                        r_transitionType = Lambda * rand(generator); //random variable for weighted sampling of transitions
                        if(r_transitionType<Beta) //S->I
                        
                            m = (int) SI * rand(generator); //transition m
                            // Add infected node to lists:
                            isInfected[si_s[m]]=true;
                            infected.push_back(si_s[m]);
                            I++;
                            Mu+=mu;
                        }
                        else //I->R
                        {
                            m = (int) I * rand(generator); //transition m
                            isInfected[infected[m]]=false;
                            // Remove drawn element from infected:
                            infected[m]=infected.back();
                            infected.pop_back();
                            I--;
                            Mu-=mu;
                        }
                        // Redo list of S->I transitions:
                        si_s.clear();
                        for(contact_iterator=(*contactList_iterator).begin(); contact_iterator!=(*contactList_iterator).end(); contact_iterator++)
                        {
                            i=(*contact_iterator).first;
                            j=(*contact_iterator).second;
                            if(isInfected[i])
                            {
                                if(!isInfected[j])
                                {
                                    si_s.push_back(j);
                                }
                            }
                            else
                            {
                                if(isInfected[j])
                                {
                                    si_s.push_back(i);
                                }
                            }
                        }
                        SI = si_s.size();
                        Beta = (double)SI*beta;
                        Lambda = Beta+Mu; //new cumulative transition rate
                        // Draw new renormalized waiting time:
                        tau = randexp(generator);
                    }
                }
                // Stop if I=0:
                if(I==0)
                {
                    stopped++;
                    break;
                }
                // read out I and R if t is divisible by res_t
                if(t % outputTimeResolution ==0)
                {
                    if(t>=T_simulation)
                    {
                        break;
                    }
                    else
                    {
                        sumI_t[q][t/outputTimeResolution] = I;
                        sumSI_t[q][t/outputTimeResolution] = SI;
                    }
                }
                t++;
            }
            t_infectionStart = 0;
        }
        hist_I[q] = I;
    }

    //-------------------------------------------------------------------------------------
    // Save epidemic data to disk:
    //-------------------------------------------------------------------------------------
    double t_simu = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

    start=std::clock();

    double t_write = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

    if (verbose)
    {
        std::cout << std::endl << "temporal Gillespie---homogeneous & Poissonian SIS: N=" << N << ", T=" << T_data << ", beta=" << beta << ", mu=" << mu << ", resolution = " << outputTimeResolution << std::endl;
        std::cout << "Simulation time: " << t_simu << ", Stopped: " << stopped << "/" << ensembleSize << std::endl;
        std::cout << "Writing to file: " << t_write << std::endl;
    }

    SI_result result;

    result.I = sumI_t;
    result.SI = sumSI_t;
    result.hist = hist_I;

    return result;
}

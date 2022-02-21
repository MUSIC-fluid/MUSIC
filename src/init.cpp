// Copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "./util.h"
#include "./cell.h"
#include "./grid.h"
#include "./init.h"
#include "./eos.h"

#ifndef _OPENMP
    #define omp_get_thread_num() 0
    #define omp_get_num_threads() 1
#else
    #include <omp.h>
#endif

using std::vector;
using std::ifstream;
using Util::hbarc;


Init::Init(const EOS &eosIn, InitData &DATA_in,
           std::shared_ptr<HydroSourceBase> hydro_source_ptr_in) :
    DATA(DATA_in), eos(eosIn){
    hydro_source_terms_ptr = hydro_source_ptr_in;
}


void Init::InitArena(SCGrid &arena_prev, SCGrid &arena_current,
                     SCGrid &arena_future) {
    print_num_of_threads();
    music_message.info("initArena");
    if (DATA.Initial_profile == 0) {
        music_message << "Using Initial_profile=" << DATA.Initial_profile;
        music_message << "nx=" << DATA.nx << ", ny=" << DATA.ny;
        music_message << "dx=" << DATA.delta_x << ", dy=" << DATA.delta_y;
        music_message.flush("info");
    } else if (DATA.Initial_profile == 1) {
        music_message << "Using Initial_profile=" << DATA.Initial_profile;
        DATA.nx = 2;
        DATA.ny = 2;
        DATA.neta = 695;
        DATA.delta_x = 0.1;
        DATA.delta_y = 0.1;
        DATA.delta_eta = 0.02;
        music_message << "nx=" << DATA.nx << ", ny=" << DATA.ny;
        music_message << "dx=" << DATA.delta_x << ", dy=" << DATA.delta_y;
        music_message << "neta=" << DATA.neta << ", deta=" << DATA.delta_eta;
        music_message.flush("info");
    } else if (DATA.Initial_profile == 8) {
        music_message.info(DATA.initName);
        ifstream profile(DATA.initName.c_str());
        if (!profile.is_open()) {
            music_message << "Initial profile: " << DATA.initName
                          << " not found.";
            music_message.flush("error");
            exit(1);
        }
        std::string dummy;
        int nx, ny, neta;
        double deta, dx, dy, dummy2;
        // read the first line with general info
        profile >> dummy >> dummy >> dummy2
                >> dummy >> neta >> dummy >> nx >> dummy >> ny
                >> dummy >> deta >> dummy >> dx >> dummy >> dy;
        profile.close();
        music_message << "Using Initial_profile=" << DATA.Initial_profile
                      << ". Overwriting lattice dimensions:";
        DATA.nx = nx;
        DATA.ny = ny;
        DATA.delta_x = dx;
        DATA.delta_y = dy;

        music_message << " neta=" << neta << ", nx=" << nx << ", ny=" << ny;
        music_message << " deta=" << DATA.delta_eta << ", dx=" << DATA.delta_x
                      << ", dy=" << DATA.delta_y;
        music_message.flush("info");
    } else if (   DATA.Initial_profile == 9 || DATA.Initial_profile == 91
               || DATA.Initial_profile == 92 || DATA.Initial_profile == 93) {
        music_message.info(DATA.initName);
        ifstream profile(DATA.initName.c_str());
        if (!profile.is_open()) {
            music_message << "Initial profile: " << DATA.initName
                          << " not found.";
            music_message.flush("error");
            exit(1);
        }
        std::string dummy;
        int nx, ny, neta;
        double deta, dx, dy, dummy2;
        // read the first line with general info
        profile >> dummy >> dummy >> dummy2
                >> dummy >> neta >> dummy >> nx >> dummy >> ny
                >> dummy >> deta >> dummy >> dx >> dummy >> dy;
        profile.close();
        music_message << "Using Initial_profile=" << DATA.Initial_profile
                      << ". Overwriting lattice dimensions:";
        DATA.nx = nx;
        DATA.ny = ny;
        DATA.neta = neta;
        DATA.delta_x = dx;
        DATA.delta_y = dy;
        DATA.delta_eta = 0.1;

        music_message << "neta=" << neta << ", nx=" << nx << ", ny=" << ny;
        music_message << "deta=" << DATA.delta_eta << ", dx=" << DATA.delta_x
                      << ", dy=" << DATA.delta_y;
        music_message.flush("info");
    } else if (DATA.Initial_profile == 94  ||
               DATA.Initial_profile == 941 ||
               DATA.Initial_profile == 95){
        //Initial profile == 94 => Bulk pressure = eps/3 - pressure => Tmunu is traceless
        //Initial profile == 941 => Bulk pressure = eps/3 - pressure => Tmunu is traceless
                              //--> Bullet is present
        //Initial profile == 95 => Bulk pressure = -pressure => Tmunu is not traceless

        //Open root file
        #ifdef ROOT_FOUND
        TFile* fIC = TFile::Open(DATA.initName.data(),"READ");
        TH3D* hE = (TH3D*) fIC->Get("energy_density");

        music_message << "Using Initial_profile=" << DATA.Initial_profile
                      << ". Overwriting lattice dimensions:";

        DATA.nx = hE->GetXaxis()->GetNbins();
        DATA.delta_x = hE->GetXaxis()->GetBinWidth(1);
        DATA.x_size = -2*(hE->GetXaxis()->GetBinCenter(1));

        DATA.ny = hE->GetYaxis()->GetNbins();
        DATA.delta_y = hE->GetYaxis()->GetBinWidth(1);
        DATA.y_size = -2*(hE->GetYaxis()->GetBinCenter(1));

        DATA.neta = hE->GetZaxis()->GetNbins();
        DATA.delta_eta = hE->GetZaxis()->GetBinWidth(1);
        DATA.eta_size = -2*(hE->GetZaxis()->GetBinCenter(1));

        music_message << " neta=" << DATA.neta << ", nx=" << DATA.nx << ", ny=" << DATA.ny;
        music_message << " deta=" << DATA.delta_eta << ", dx=" << DATA.delta_x
                      << ", dy=" << DATA.delta_y;
        music_message << "Leta= "<< DATA.eta_size <<", Lx=" << DATA.x_size<< ", Ly="<<DATA.y_size;
        music_message.flush("info");

        fIC->Close();
        #else
        music_message << "MUSIC was not compiled against ROOT, necessary for ";
        music_message << "reading ICs with Initial profile = 94, 941 or 95.";
        music_message.flush("error");
        #endif
    } else if (DATA.Initial_profile == 11 || DATA.Initial_profile == 111) {
        double tau_overlap = 2.*7./(sinh(DATA.beam_rapidity));
        DATA.tau0 = std::max(DATA.tau0, tau_overlap);
        music_message << "tau0 = " << DATA.tau0 << " fm/c.";
        music_message.flush("info");
    } else if (DATA.Initial_profile == 112) {
        double tau_overlap = 2.*7./(sinh(DATA.beam_rapidity));
        DATA.tau0 = std::max(DATA.tau0, tau_overlap) - DATA.delta_tau;
        music_message << "tau0 = " << DATA.tau0 << " fm/c.";
        music_message.flush("info");
    } else if (DATA.Initial_profile == 13 || DATA.Initial_profile == 131) {
        DATA.tau0 = (hydro_source_terms_ptr.lock()->get_source_tau_min()
                     - DATA.delta_tau);
        DATA.tau0 = static_cast<int>(DATA.tau0/0.02)*0.02;
        DATA.tau0 = std::max(0.1, DATA.tau0);
    } else if (DATA.Initial_profile == 30) {
        DATA.tau0 = hydro_source_terms_ptr.lock()->get_source_tau_min();
    } else if (DATA.Initial_profile == 42) {
        // initial condition from the JETSCAPE framework
        music_message << "Using Initial_profile=" << DATA.Initial_profile
                      << ". Overwriting lattice dimensions:";
        music_message.flush("info");

        const int nx = static_cast<int>(
                sqrt(jetscape_initial_energy_density.size()/DATA.neta));
        const int ny = nx;
        DATA.nx = nx;
        DATA.ny = ny;
        DATA.x_size = DATA.delta_x*nx;
        DATA.y_size = DATA.delta_y*ny;

        music_message << "neta = " << DATA.neta
                      << ", nx = " << nx << ", ny = " << ny;
        music_message.flush("info");
        music_message << "deta=" << DATA.delta_eta
                      << ", dx=" << DATA.delta_x
                      << ", dy=" << DATA.delta_y;
        music_message.flush("info");
        music_message << "x_size = "     << DATA.x_size
                      << ", y_size = "   << DATA.y_size
                      << ", eta_size = " << DATA.eta_size;
        music_message.flush("info");
    } else if (DATA.Initial_profile == 101) {
        music_message << "Using Initial_profile = " << DATA.Initial_profile;
        music_message.flush("info");
        music_message << "nx = " << DATA.nx << ", ny = " << DATA.ny
                      << ", neta = " << DATA.neta;
        music_message.flush("info");
        music_message << "dx = " << DATA.delta_x << ", dy = " << DATA.delta_y
                      << ", deta = " << DATA.delta_eta;
        music_message.flush("info");
    } else if (DATA.Initial_profile == 17) {
        music_message << "Using NeXus initial profile. (Initial_profile = "
                      << DATA.Initial_profile <<")";
        music_message.flush("info");
        music_message << "Reading IC from "<< DATA.initName.c_str();
        music_message.flush("info");
        //Open nexus file
        ifstream profile_nexus(DATA.initName.c_str());
        //Read header
        if (read_nexus_header(&profile_nexus) < 0 ) return;
        profile_nexus.close();
        music_message << "Overwriting lattice dimensions: ";
        music_message << "neta = " << DATA.neta
                      << ", nx = " << DATA.nx << ", ny = " << DATA.ny;
        music_message.flush("info");
        music_message << "deta=" << DATA.delta_eta
                      << ", dx=" << DATA.delta_x
                      << ", dy=" << DATA.delta_y;
        music_message.flush("info");
        music_message << "x_size = "     << DATA.x_size
                      << ", y_size = "   << DATA.y_size
                      << ", eta_size = " << DATA.eta_size;
        music_message.flush("info");

    } else if (DATA.Initial_profile == 22){
        AMPT_smeared_header();
    }

    // initialize arena
    arena_prev    = SCGrid(DATA.nx, DATA.ny, DATA.neta);
    arena_current = SCGrid(DATA.nx, DATA.ny, DATA.neta);
    arena_future  = SCGrid(DATA.nx, DATA.ny, DATA.neta);
    music_message.info("Grid allocated.");

    InitTJb(arena_prev, arena_current);

    if (DATA.output_initial_density_profiles == 1) {
        output_initial_density_profiles(arena_current);
    }
}/* InitArena */


void Init::print_num_of_threads() {
    #pragma omp parallel for
    for (int i = 0; i < 2; i++) {
        if (i == 0) {
            music_message << "OpenMP: using " << omp_get_num_threads()
                          << " threads.";
            music_message.flush("info");
        }
    }
}
void Init::AMPT_smeared_header(){

    ifstream input_file(DATA.initName.c_str());
    std::string input_line;

    music_message << "Overwritting grid dimensions..."
                  << "\n----------\n";

    auto return_val = [](std::string line, std::string par){
            size_t found = line.find(par);
            if (found!=std::string::npos){
                return line.substr(found+par.size(),line.size() );
            } else {
                return std::string("");
            }
    };

    while (getline(input_file, input_line)) {
        if (input_line[0] == '#') { // read only header lines
            if (return_val(input_line,"nx =").compare("")){
                DATA.nx = std::stoi(return_val(input_line,"nx ="));
                music_message << "nx = " << DATA.nx << "\n";
            }
            if (return_val(input_line,"ny =").compare("")){
                DATA.ny = std::stoi(return_val(input_line,"ny ="));
                music_message << "ny = " << DATA.ny << "\n";
            }
            if (return_val(input_line,"neta =").compare("")){
                DATA.neta = std::stoi(return_val(input_line,"neta ="));
                music_message << "neta = " << DATA.neta << "\n";
            }
            if (return_val(input_line,"Lx =").compare("")){
                DATA.x_size = std::stod(return_val(input_line,"Lx ="));
                music_message << "Lx = " << DATA.x_size << "\n";
            }
            if (return_val(input_line,"Ly =").compare("")){
                DATA.y_size = std::stod(return_val(input_line,"Ly ="));
                music_message << "Ly = " << DATA.y_size << "\n";
            }
            if (return_val(input_line,"Leta =").compare("")){
                DATA.eta_size = std::stod(return_val(input_line,"Leta ="));
                music_message << "Leta = " << DATA.eta_size << "\n";
            }
        }
    }

    //Now we determine grid spacing
    DATA.delta_x = DATA.x_size/(DATA.nx-1);
    DATA.delta_y = DATA.y_size/(DATA.ny-1);
    DATA.delta_eta = DATA.eta_size/(DATA.neta-1);
    music_message <<   "dx = "<<DATA.delta_x
                  << "\ndy = "<<DATA.delta_y
                  << "\ndeta = "<<DATA.delta_eta;


    music_message << "\n----------";
    music_message.flush("info");

    input_file.close();


}


///////////////////////////////////////////////////////////////////////////////
///\author Willian M. Serenone
///\brief Reads the initial condition from smeared AMPT event. Format should be
/// x y eta epsilon ux uy ueta trace pitautau pitaux pitauy pitaueta pixx pixy pixeta piyy piyeta pietaeta rhob
void Init::initial_AMPT_smeared(SCGrid &arena_prev, SCGrid &arena_current){

    const size_t nx = arena_current.nX();
    const size_t ny = arena_current.nY();
    const size_t neta = arena_current.nEta();

    const float dx = DATA.delta_x;
    const float dy = DATA.delta_y;
    const float deta = DATA.delta_eta;

    const float Lx = DATA.x_size;
    const float Ly = DATA.y_size;
    const float Leta = DATA.eta_size;

    ifstream input_file(DATA.initName.c_str());
    std::string input_line;

    long ii=0;
    while (getline(input_file, input_line)) {
        if (input_line[0] != '#') { // ignore header lines
            std::istringstream line_stream(input_line);
            std::vector<double> line_vector((std::istream_iterator<double>(line_stream)), std::istream_iterator<double>());

            //Check if the number of columns is the one expected
            if (line_vector.size() != 19){
                music_message << "Wrong number of columns in input file";
                music_message.flush("error");
                exit(1);
            }

            double x        = line_vector[0];
            double y        = line_vector[1];
            double eta      = line_vector[2];
            double eps      = line_vector[3]/Util::hbarc;
            double ux       = line_vector[4];
            double uy       = line_vector[5];
            double ueta     = line_vector[6];
            double trace    = line_vector[7]/Util::hbarc;
            double pitautau = line_vector[8]/Util::hbarc;
            double pitaux   = line_vector[9]/Util::hbarc;
            double pitauy   = line_vector[10]/Util::hbarc;
            double pitaueta = line_vector[11]/Util::hbarc;
            double pixx     = line_vector[12]/Util::hbarc;
            double pixy     = line_vector[13]/Util::hbarc;
            double pixeta   = line_vector[14]/Util::hbarc;
            double piyy     = line_vector[15]/Util::hbarc;
            double piyeta   = line_vector[16]/Util::hbarc;
            double pietaeta = line_vector[17]/Util::hbarc;
            double rhob     = line_vector[18];

            int ix = std::round( (x + Lx*.5)/dx);
            int iy = std::round( (y + Ly*.5)/dy);
            int ieta = std::round( (eta + Leta*.5)/deta);

            if((ix < 0) || (ix >= DATA.nx)){
                music_message << "ix = " << ix << " is out of boundaries";
                music_message.flush("error");
                music_message << "x = " << x;
                music_message.flush("error");
                exit(1);
            }
            if((iy < 0) || (iy >= DATA.ny)){
                music_message << "iy = " << iy << " is out of boundaries";
                music_message.flush("error");
                music_message << "y = " << y;
                music_message.flush("error");
                exit(1);
            }
            if((ieta < 0) || (ieta >= DATA.neta)){
                music_message << "ieta = " << ieta << " is out of boundaries";
                music_message.flush("error");
                music_message << "eta_s = " << eta;
                music_message.flush("error");
                exit(1);
            }


            double tau_ueta = DATA.tau0*ueta;
            double utau = sqrt(1+ux*ux + uy*uy + tau_ueta*tau_ueta);

            arena_current(ix, iy, ieta).epsilon = eps;
            arena_current(ix, iy, ieta).rhob = rhob;
            arena_current(ix, iy, ieta).u[0] = utau;
            arena_current(ix, iy, ieta).u[1] = ux;
            arena_current(ix, iy, ieta).u[2] = uy;
            arena_current(ix, iy, ieta).u[3] = tau_ueta;

            arena_current(ix, iy, ieta).Wmunu[0] = pitautau;
            arena_current(ix, iy, ieta).Wmunu[1] = pitaux;
            arena_current(ix, iy, ieta).Wmunu[2] = pitauy;
            arena_current(ix, iy, ieta).Wmunu[3] = DATA.tau0*pitaueta;
            arena_current(ix, iy, ieta).Wmunu[4] = pixx;
            arena_current(ix, iy, ieta).Wmunu[5] = pixy;
            arena_current(ix, iy, ieta).Wmunu[6] = DATA.tau0*pixeta;
            arena_current(ix, iy, ieta).Wmunu[7] = piyy;
            arena_current(ix, iy, ieta).Wmunu[8] = DATA.tau0*piyeta;
            arena_current(ix, iy, ieta).Wmunu[9] = DATA.tau0*DATA.tau0*pietaeta;

            arena_current(ix, iy, ieta).pi_b = (eps-trace)/3.0 - eos.get_pressure(eps,rhob);

            arena_prev(ix, iy, ieta) = arena_current(ix, iy, ieta);
        }
    }
    input_file.close();
}


void Init::initial_trento_XY(int ieta, SCGrid &arena_prev, SCGrid &arena_current) {
    // initial condition is a 2D profile generated by Trento
    const size_t nx = arena_current.nX();
    const size_t ny = arena_current.nY();
    ifstream surfaceFile(DATA.initName.c_str());
    size_t iy = 0;
    std::string surfaceLine;
    int entropy_flag = DATA.initializeEntropy;

    //Variables for eta envelope
    double eta = (DATA.delta_eta)*ieta - (DATA.eta_size)/2.0;
    double eta_envelop_ed = eta_profile_plateau(eta, DATA.eta_flat/2.0,
                                                DATA.eta_fall_off);

    while (getline(surfaceFile, surfaceLine)) {
        if (surfaceLine[0] != '#') { // ignore TRENTo header lines
            std::istringstream lineStream(surfaceLine);
            vector<double> lineVector((std::istream_iterator<double>(lineStream)), std::istream_iterator<double>());
            // checks if number of columns in initial condition file matches grid nx size
            if (lineVector.size() != nx) {
                music_message.error("nx size on initial condition file does not match MUSIC nx !");
                exit(1);
            }
            for (size_t ix = 0; ix < nx; ix++) {
               double s = lineVector[ix];
               double epsilon;
               if (entropy_flag == 1) {
                 epsilon = eos.get_s2e(s, 0.0);
               } else if (entropy_flag == 0) {
                 epsilon = s;
               }
               epsilon = std::max(epsilon, Util::small_eps);
               double rhob = 0.;
               // set all values in the grid element:
               arena_current(ix, iy, ieta).epsilon = epsilon*eta_envelop_ed*DATA.sFactor;
               arena_current(ix, iy, ieta).rhob = rhob;
               /* for HIC */
                arena_current(ix, iy, ieta).u[0] = 1.0;
                arena_current(ix, iy, ieta).u[1] = 0.0;
                arena_current(ix, iy, ieta).u[2] = 0.0;
                arena_current(ix, iy, ieta).u[3] = 0.0;

                //Shear tensor is zero
                arena_current(ix, iy, ieta).Wmunu[0] = 0.0;
                arena_current(ix, iy, ieta).Wmunu[1] = 0.0;
                arena_current(ix, iy, ieta).Wmunu[2] = 0.0;
                arena_current(ix, iy, ieta).Wmunu[3] = 0.0;
                arena_current(ix, iy, ieta).Wmunu[4] = 0.0;
                arena_current(ix, iy, ieta).Wmunu[5] = 0.0;
                arena_current(ix, iy, ieta).Wmunu[6] = 0.0;
                arena_current(ix, iy, ieta).Wmunu[7] = 0.0;
                arena_current(ix, iy, ieta).Wmunu[8] = 0.0;
                arena_current(ix, iy, ieta).Wmunu[9] = 0.0;


                //Init bulk
                if (DATA.turn_on_bulk == 0){
                    arena_current(ix, iy, ieta).pi_b = 0.0;
                } else {
                    if (DATA.Initial_profile == 12){
                        arena_current(ix, iy, ieta).pi_b = 0.0;
                    } else if (DATA.Initial_profile == 121){
                        double pressure = eos.get_pressure(epsilon, rhob);
                        arena_current(ix, iy, ieta).pi_b = epsilon/3. - pressure;
                    }
                }




                arena_prev(ix, iy, ieta) = arena_current(ix, iy, ieta);
            } /* x */
            iy++;
        }
    } // end of initial condition file
    if (iy != ny) {
        music_message.error("ny size on initial condition file does not match MUSIC nx!");
        exit(1);
    }
}

int Init::read_nexus_header(ifstream* nexus_IC){
    std::string dummy;
    *nexus_IC >> dummy; //energy
    *nexus_IC >> dummy >> dummy >> dummy >> dummy; //zproj mproj ztarg mtarg
    *nexus_IC >> dummy; //bmax
    *nexus_IC >> dummy >> dummy >> dummy;
    *nexus_IC >> DATA.nx >> DATA.ny >> DATA.neta;
    double xmin, xmax, ymin, ymax, eta_min, eta_max;
    *nexus_IC >> xmin >> xmax >> ymin >> ymax >> eta_min >> eta_max;

    if(xmin > xmax){
        music_message << "xmin > xmax. Check NeXus file format.";
        music_message << "xmin = " << xmin<<". ";
        music_message << "xmax = " << xmax<<".";
        music_message.flush("error");
        return - 1;
    }

    if(ymin > ymax){
        music_message << "ymin > ymax. Check NeXus file format.";
        music_message << "ymin = " << ymin<<". ";
        music_message << "ymax = " << ymax<<".";
        music_message.flush("error");
        return -2;
    }

    if(eta_min > eta_max){
        music_message << "eta_min > eta_max. Check NeXus file format.";
        music_message.flush("error");
        music_message << "eta_min = " << eta_min<<". ";
        music_message << "eta_max = " << eta_max<<".";
        music_message.flush("error");
        return -3;
    }

    if(abs(xmin+xmax) > 1.E-10){
        music_message << "Grid size in x is asymmetrical. We will symmetrize it";
        music_message.flush("warning");
    }
    if(abs(ymin+ymax) > 1.E-10){
        music_message << "Grid size in x is asymmetrical. We will symmetrize it";
        music_message.flush("warning");
    }
    if(abs(eta_min+eta_max) > 1.E-10){
        music_message << "Grid size in eta_s is asymmetrical. We will symmetrize it";
        music_message.flush("warning");
    }

    DATA.x_size = xmax-xmin;
    DATA.y_size = ymax-ymin;
    DATA.eta_size = eta_max-eta_min;

    DATA.delta_x = DATA.x_size/DATA.nx;
    DATA.delta_y = DATA.y_size/DATA.ny;
    DATA.delta_eta = DATA.eta_size/DATA.neta;

    return 0;

}

int Init::read_nexus_profile(SCGrid &arena_prev,
                             SCGrid &arena_current){

    ifstream profile_nexus(DATA.initName.c_str());
    if (read_nexus_header(&profile_nexus) < 0 ) return -1;
    double tau = DATA.tau0;
    std::string dummy;


    const int nx = arena_current.nX();
    const int ny = arena_current.nY();
    const int neta = arena_current.nEta();

    std::vector<std::vector<std::vector<double>>> temp_profile_ed =
        std::vector<std::vector<std::vector<double>>>(neta,std::vector<std::vector<double>>(ny,std::vector<double>(nx,.0)));
    std::vector<std::vector<std::vector<double>>> temp_profile_vx =
        std::vector<std::vector<std::vector<double>>>(neta,std::vector<std::vector<double>>(ny,std::vector<double>(nx,.0)));
    std::vector<std::vector<std::vector<double>>> temp_profile_vy =
        std::vector<std::vector<std::vector<double>>>(neta,std::vector<std::vector<double>>(ny,std::vector<double>(nx,.0)));
    std::vector<std::vector<std::vector<double>>> temp_profile_veta =
        std::vector<std::vector<std::vector<double>>>(neta,std::vector<std::vector<double>>(ny,std::vector<double>(nx,.0)));

    std::vector<std::vector<std::vector<double>>> temp_profile_net_up =
        std::vector<std::vector<std::vector<double>>>(neta,std::vector<std::vector<double>>(ny,std::vector<double>(nx,.0)));
    std::vector<std::vector<std::vector<double>>> temp_profile_net_down =
        std::vector<std::vector<std::vector<double>>>(neta,std::vector<std::vector<double>>(ny,std::vector<double>(nx,.0)));
    std::vector<std::vector<std::vector<double>>> temp_profile_net_strange =
        std::vector<std::vector<std::vector<double>>>(neta,std::vector<std::vector<double>>(ny,std::vector<double>(nx,.0)));

    // Ignore T^\mu \nu
    for (int ignore=0; ignore<16*nx*ny*neta;++ignore) profile_nexus >> dummy;

    // Ignore flavor current
    for (int ignore=0; ignore<12*nx*ny*neta;++ignore) profile_nexus >> dummy;

    // Grab energy density
    for (int ieta = 0; ieta < neta; ++ieta)
        for (int iy = 0; iy < ny; ++iy)
            for (int ix = 0; ix < nx; ++ix)
                profile_nexus >> temp_profile_ed[ieta][iy][ix];

    // Grab velocity (non-relativistic)
    for (int ieta = 0; ieta < neta; ++ieta)
        for (int iy = 0; iy < ny; ++iy)
            for (int ix = 0; ix < nx; ++ix)
                profile_nexus >> temp_profile_vx[ieta][iy][ix]
                              >> temp_profile_vy[ieta][iy][ix]
                              >> temp_profile_veta[ieta][iy][ix];

    // Grab net flavor densities
    for (int ieta = 0; ieta < neta; ++ieta)
        for (int iy = 0; iy < ny; ++iy)
            for (int ix = 0; ix < nx; ++ix){
                profile_nexus >> temp_profile_net_up[ieta][iy][ix];
                profile_nexus >> temp_profile_net_down[ieta][iy][ix];
                profile_nexus >> temp_profile_net_strange[ieta][iy][ix];
            }

    for (int ieta = 0; ieta < neta; ++ieta){
        for (int iy = 0; iy< ny; iy++) {
            for (int ix = 0; ix < nx; ix++) {

                double rhob = (temp_profile_net_up[ieta][iy][ix]
                               + temp_profile_net_down[ieta][iy][ix]
                               + temp_profile_net_strange[ieta][iy][ix])/3.;
                double epsilon = (temp_profile_ed[ieta][iy][ix])/hbarc; // Convert to 1/fm^4
                if (epsilon < 1.E-12) epsilon = 1.E-12;

                if ( isnan(epsilon) ){
                    music_message << "nan found in energy profile";
                    music_message.flush("error");
                }
                if ( isnan(temp_profile_vx[ieta][iy][ix]) ){
                    music_message << "nan found in vx profile";
                    music_message.flush("error");
                }
                if ( isnan(temp_profile_vy[ieta][iy][ix]) ){
                    music_message << "nan found in vy profile";
                    music_message.flush("error");
                }
                if ( isnan(temp_profile_veta[ieta][iy][ix]) ){
                    music_message << "nan found in veta profile";
                    music_message.flush("error");
                }


                arena_current(ix, iy, ieta).epsilon = epsilon;
                arena_current(ix, iy, ieta).rhob = rhob;

                double vx = temp_profile_vx[ieta][iy][ix];
                double vy = temp_profile_vy[ieta][iy][ix];
                double veta = temp_profile_veta[ieta][iy][ix];

                double v2 = vx*vx + vy*vy + pow(tau*veta,2);
                if (v2 > 1){
                    music_message << "Velocity greater than light detected in IC";
                    music_message.flush("error");
                }

                double gamma = pow(1.- v2,.5);


                arena_current(ix, iy, ieta).u[0] = gamma;
                arena_current(ix, iy, ieta).u[1] = gamma*vx;
                arena_current(ix, iy, ieta).u[2] = gamma*vy;
                arena_current(ix, iy, ieta).u[3] = gamma*veta;



                if (DATA.turn_on_bulk == 1) {
                    //Make the tensor traceless
                    double pressure = eos.get_pressure(epsilon, rhob);
                    arena_current(ix, iy, ieta).pi_b = epsilon/3. - pressure;
                } else {
                    arena_current(ix, iy, ieta).pi_b = 0;
                }


                //Setup viscous part
                double temp_profile_pixx   = 0.0;
                double temp_profile_pixy   = 0.0;
                double temp_profile_pixeta = 0.0;
                double temp_profile_piyy   = 0.0;
                double temp_profile_piyeta = 0.0;

                double utau = gamma;
                double ueta = tau*gamma*veta;
                double ux = vx*gamma;
                double uy = vy*gamma;
                double temp_profile_pietaeta = (
                (2.*(  ux*uy*temp_profile_pixy
                     + ux*ueta*temp_profile_pixeta
                     + uy*ueta*temp_profile_piyeta)
                 - (utau*utau - ux*ux)*temp_profile_pixx
                 - (utau*utau - uy*uy)*temp_profile_piyy)
                /(utau*utau - ueta*ueta));
                double temp_profile_pitaux   = (1./utau
                   *(  temp_profile_pixx*ux + temp_profile_pixy*uy
                     + temp_profile_pixeta*ueta));
                double temp_profile_pitauy  = (1./utau
                *(  temp_profile_pixy*ux
                  + temp_profile_piyy*uy
                  + temp_profile_piyeta*ueta));
                double temp_profile_pitaueta = (1./utau
                *(  temp_profile_pixeta*ux
                  + temp_profile_piyeta*uy
                  + temp_profile_pietaeta*ueta));
                double temp_profile_pitautau = (1./utau
                *(  temp_profile_pitaux*ux
                  + temp_profile_pitauy*uy
                  + temp_profile_pitaueta*ueta));

                arena_current(ix, iy, ieta).Wmunu[0] = temp_profile_pitautau;
                arena_current(ix, iy, ieta).Wmunu[1] = temp_profile_pitaux;
                arena_current(ix, iy, ieta).Wmunu[2] = temp_profile_pitauy;
                arena_current(ix, iy, ieta).Wmunu[3] = temp_profile_pitaueta;
                arena_current(ix, iy, ieta).Wmunu[4] = temp_profile_pixx;
                arena_current(ix, iy, ieta).Wmunu[5] = temp_profile_pixy;
                arena_current(ix, iy, ieta).Wmunu[6] = temp_profile_pixeta;
                arena_current(ix, iy, ieta).Wmunu[7] = temp_profile_piyy;
                arena_current(ix, iy, ieta).Wmunu[8] = temp_profile_piyeta;
                arena_current(ix, iy, ieta).Wmunu[9] = temp_profile_pietaeta;

                arena_prev(ix, iy, ieta) = arena_current(ix, iy, ieta);
            }
        }
    }
    profile_nexus.close();
    return 0;
}

//! This is a shell function to initial hydrodynamic fields
void Init::InitTJb(SCGrid &arena_prev, SCGrid &arena_current) {
    if (DATA.Initial_profile == 0) {
        // Gubser flow test
        music_message.info(" Perform Gubser flow test ... ");
        music_message.info(" ----- information on initial distribution -----");

        #pragma omp parallel for
        for (int ieta = 0; ieta < arena_current.nEta(); ieta++) {
            initial_Gubser_XY(ieta, arena_prev, arena_current);
        }
    } else if (DATA.Initial_profile == 1) {
        // code test in 1+1 D vs Monnai's results
        music_message.info(" Perform 1+1D test vs Monnai's results... ");
        initial_1p1D_eta(arena_prev, arena_current);
    } else if (DATA.Initial_profile == 8) {
        // read in the profile from file
        // - IPGlasma initial conditions with initial flow
        music_message.info(" ----- information on initial distribution -----");
        music_message << "file name used: " << DATA.initName;
        music_message.flush("info");

        #pragma omp parallel for
        for (int ieta = 0; ieta < arena_current.nEta(); ieta++) {
            initial_IPGlasma_XY(ieta, arena_prev, arena_current);
        }
    } else if (   DATA.Initial_profile == 9 || DATA.Initial_profile == 91
               || DATA.Initial_profile == 92 || DATA.Initial_profile == 93) {
        // read in the profile from file
        // - IPGlasma initial conditions with initial flow
        // and initial shear viscous tensor
        music_message.info(" ----- information on initial distribution -----");
        music_message << "file name used: " << DATA.initName;
        music_message.flush("info");

        #pragma omp parallel for
        for (int ieta = 0; ieta < arena_current.nEta(); ieta++) {
            initial_IPGlasma_XY_with_pi(ieta, arena_prev, arena_current);
        }
    } else if (   DATA.Initial_profile == 94  ||
                  DATA.Initial_profile == 941 ||
                  DATA.Initial_profile == 95) {
        // read in the profile from file
        // - Full Tmunu in ROOT format
        music_message.info(" ----- information on initial distribution -----");
        music_message << "file name used: " << DATA.initName;
        music_message.flush("info");

        initial_ROOT_XYEta_with_pi(arena_prev, arena_current);

    } else if (DATA.Initial_profile == 11 || DATA.Initial_profile == 111) {
        // read in the transverse profile from file with finite rho_B
        // the initial entropy and net baryon density profile are
        // constructed by nuclear thickness function TA and TB.
        // Along the longitudinal direction an asymmetric contribution from
        // target and projectile thickness function is allowed
        music_message.info(" ----- information on initial distribution -----");
        music_message << "file name used: " << DATA.initName_TA << " and "
                      << DATA.initName_TB;
        music_message.flush("info");

        initial_MCGlb_with_rhob(arena_prev, arena_current);
    } else if (DATA.Initial_profile == 112) {
        music_message.info(
                "Initialize hydro with source terms from TA and TB");
        #pragma omp parallel for
        for (int ieta = 0; ieta < arena_current.nEta(); ieta++) {
            initial_with_zero_XY(ieta, arena_prev, arena_current);
        }
    } else if (DATA.Initial_profile == 13 || DATA.Initial_profile == 131) {
        music_message.info("Initialize hydro with source terms");
        #pragma omp parallel for
        for (int ieta = 0; ieta < arena_current.nEta(); ieta++) {
            initial_with_zero_XY(ieta, arena_prev, arena_current);
        }
    } else if (DATA.Initial_profile == 30) {
        #pragma omp parallel for
        for (int ieta = 0; ieta < arena_current.nEta(); ieta++) {
            initial_AMPT_XY(ieta, arena_prev, arena_current);
        }
    } else if (DATA.Initial_profile == 42) {
        // initialize hydro with vectors from JETSCAPE
        music_message.info(" ----- information on initial distribution -----");
        music_message << "initialized with a JETSCAPE initial condition.";
        music_message.flush("info");
        #pragma omp parallel for
        for (int ieta = 0; ieta < arena_current.nEta(); ieta++) {
            initial_with_jetscape(ieta, arena_prev, arena_current);
        }
        clean_up_jetscape_arrays();
    } else if (DATA.Initial_profile == 101) {
        music_message.info(" ----- information on initial distribution -----");
        music_message << "file name used: " << DATA.initName;
        music_message.flush("info");
        initial_UMN_with_rhob(arena_prev, arena_current);
    } else if (DATA.Initial_profile == 17) {
        music_message.info(" ----- information on initial distribution -----");
        music_message << "file name used: " << DATA.initName;
        music_message.flush("info");
        read_nexus_profile(arena_prev, arena_current);
        music_message << "Completed " << DATA.initName;
        music_message.flush("info");
    } else if (DATA.Initial_profile == 12 || DATA.Initial_profile == 121) {
        // reads transverse profile in 2-D TRENTo output format
        // see: http://qcd.phy.duke.edu/trento/index.html
        // 12 -> Pi_b = 0. Tmunu is not tracceless
        // 121 -> Pi_b = eps/3 - P. Tmunu is traceless
        music_message.info("Reading initial profile in TRENTo format");
        #pragma omp parallel for
        for (int ieta = 0; ieta < arena_current.nEta(); ieta++){
            initial_trento_XY(ieta, arena_prev, arena_current);
        }
    } else if (DATA.Initial_profile == 22){
        music_message.info("Reading initial profile in AMPT-genesis format");
        initial_AMPT_smeared(arena_prev, arena_current);
    }

    if (DATA.viscosity_flag == 0) {
        // for ideal hydrodynamic simualtions set all viscous tensor to zero
        music_message << "Running ideal hydrodynamic simulations ...";
        music_message.flush("info");
        music_message << "Setting the initial viscous tensor to zero.";
        music_message.flush("info");
        const int grid_neta = arena_current.nEta();
        const int grid_nx   = arena_current.nX();
        const int grid_ny   = arena_current.nY();
        #pragma omp parallel for collapse(3)
        for (int ieta = 0; ieta < grid_neta; ieta++) {
            for (int ix = 0; ix < grid_nx; ix++) {
                for (int iy = 0; iy < grid_ny; iy++) {
                    arena_prev(ix, iy, ieta).Wmunu = {0.};
                    arena_prev(ix, iy, ieta).pi_b = 0.;
                    arena_current(ix, iy, ieta).Wmunu = {0.};
                    arena_current(ix, iy, ieta).pi_b = 0.;
                }
            }
        }
    }
    music_message.info("initial distribution done.");
}

void Init::initial_Gubser_XY(int ieta, SCGrid &arena_prev,
                             SCGrid &arena_current) {
    std::string input_filename;
    std::string input_filename_prev;
    if (DATA.turn_on_shear == 1) {
        input_filename = "tests/Gubser_flow/Initial_Profile.dat";
    } else {
        input_filename = "tests/Gubser_flow/y=0_tau=1.00_ideal.dat";
        input_filename_prev = "tests/Gubser_flow/y=0_tau=0.98_ideal.dat";
    }

    ifstream profile(input_filename.c_str());
    if (!profile.good()) {
        music_message << "Init::InitTJb: "
                      << "Can not open the initial file: " << input_filename;
        music_message.flush("error");
        exit(1);
    }
    ifstream profile_prev;
    if (DATA.turn_on_shear == 0) {
        profile_prev.open(input_filename_prev.c_str());
        if (!profile_prev.good()) {
            music_message << "Init::InitTJb: "
                          << "Can not open the initial file: "
                          << input_filename_prev;
            music_message.flush("error");
            exit(1);
        }
    }

    const int nx = arena_current.nX();
    const int ny = arena_current.nY();
    double temp_profile_ed[nx][ny];
    double temp_profile_ux[nx][ny];
    double temp_profile_uy[nx][ny];
    double temp_profile_ed_prev[nx][ny];
    double temp_profile_rhob[nx][ny];
    double temp_profile_rhob_prev[nx][ny];
    double temp_profile_ux_prev[nx][ny];
    double temp_profile_uy_prev[nx][ny];
    double temp_profile_pixx[nx][ny];
    double temp_profile_piyy[nx][ny];
    double temp_profile_pixy[nx][ny];
    double temp_profile_pi00[nx][ny];
    double temp_profile_pi0x[nx][ny];
    double temp_profile_pi0y[nx][ny];
    double temp_profile_pi33[nx][ny];

    double dummy;
    for (int ix = 0; ix < nx; ix++) {
        for (int iy = 0; iy < ny; iy++) {
            if (DATA.turn_on_shear == 1) {
                profile >> dummy >> dummy >> temp_profile_ed[ix][iy]
                        >> temp_profile_ux[ix][iy] >> temp_profile_uy[ix][iy];
                profile >> temp_profile_pixx[ix][iy]
                        >> temp_profile_piyy[ix][iy]
                        >> temp_profile_pixy[ix][iy]
                        >> temp_profile_pi00[ix][iy]
                        >> temp_profile_pi0x[ix][iy]
                        >> temp_profile_pi0y[ix][iy]
                        >> temp_profile_pi33[ix][iy];
            } else {
                profile >> dummy >> dummy >> temp_profile_ed[ix][iy]
                        >> temp_profile_rhob[ix][iy]
                        >> temp_profile_ux[ix][iy] >> temp_profile_uy[ix][iy];
                profile_prev >> dummy >> dummy >> temp_profile_ed_prev[ix][iy]
                             >> temp_profile_rhob_prev[ix][iy]
                             >> temp_profile_ux_prev[ix][iy]
                             >> temp_profile_uy_prev[ix][iy];
            }
        }
    }
    profile.close();
    if (DATA.turn_on_shear == 0) {
        profile_prev.close();
    }

    for (int ix = 0; ix < nx; ix++) {
        for (int iy = 0; iy< ny; iy++) {
            double rhob = 0.0;
            if (DATA.turn_on_shear == 0 && DATA.turn_on_rhob == 1) {
                rhob = temp_profile_rhob[ix][iy];
            }

            double epsilon = temp_profile_ed[ix][iy];

            arena_current(ix, iy, ieta).epsilon = epsilon;
            arena_prev   (ix, iy, ieta).epsilon = epsilon;
            arena_current(ix, iy, ieta).rhob    = rhob;
            arena_prev   (ix, iy, ieta).rhob    = rhob;

            double utau_local = sqrt(1.
                          + temp_profile_ux[ix][iy]*temp_profile_ux[ix][iy]
                          + temp_profile_uy[ix][iy]*temp_profile_uy[ix][iy]);
            arena_current(ix, iy, ieta).u[0] = utau_local;
            arena_current(ix, iy, ieta).u[1] = temp_profile_ux[ix][iy];
            arena_current(ix, iy, ieta).u[2] = temp_profile_uy[ix][iy];
            arena_current(ix, iy, ieta).u[3] = 0.0;
            arena_prev(ix, iy, ieta).u = arena_current(ix, iy, ieta).u;

            if (DATA.turn_on_shear == 0) {
                double utau_prev = sqrt(1.
                    + temp_profile_ux_prev[ix][iy]*temp_profile_ux_prev[ix][iy]
                    + temp_profile_uy_prev[ix][iy]*temp_profile_uy_prev[ix][iy]
                );
                arena_prev(ix, iy, ieta).u[0] = utau_prev;
                arena_prev(ix, iy, ieta).u[1] = temp_profile_ux_prev[ix][iy];
                arena_prev(ix, iy, ieta).u[2] = temp_profile_uy_prev[ix][iy];
                arena_prev(ix, iy, ieta).u[3] = 0.0;
            }

            if (DATA.turn_on_shear == 1) {
                arena_current(ix,iy,ieta).Wmunu[0] = temp_profile_pi00[ix][iy];
                arena_current(ix,iy,ieta).Wmunu[1] = temp_profile_pi0x[ix][iy];
                arena_current(ix,iy,ieta).Wmunu[2] = temp_profile_pi0y[ix][iy];
                arena_current(ix,iy,ieta).Wmunu[3] = 0.0;
                arena_current(ix,iy,ieta).Wmunu[4] = temp_profile_pixx[ix][iy];
                arena_current(ix,iy,ieta).Wmunu[5] = temp_profile_pixy[ix][iy];
                arena_current(ix,iy,ieta).Wmunu[6] = 0.0;
                arena_current(ix,iy,ieta).Wmunu[7] = temp_profile_piyy[ix][iy];
                arena_current(ix,iy,ieta).Wmunu[8] = 0.0;
                arena_current(ix,iy,ieta).Wmunu[9] = temp_profile_pi33[ix][iy];
            }
            arena_prev(ix,iy,ieta).Wmunu = arena_current(ix,iy,ieta).Wmunu;
        }
    }
}

void Init::initial_1p1D_eta(SCGrid &arena_prev, SCGrid &arena_current) {
    std::string input_ed_filename;
    std::string input_rhob_filename;
    input_ed_filename = "tests/test_1+1D_with_Akihiko/e_baryon_init.dat";
    input_rhob_filename = "tests/test_1+1D_with_Akihiko/rhoB_baryon_init.dat";

    ifstream profile_ed(input_ed_filename.c_str());
    if (!profile_ed.good()) {
        music_message << "Init::InitTJb: "
                      << "Can not open the initial file: "
                      << input_ed_filename;
        music_message.flush("error");
        exit(1);
    }
    ifstream profile_rhob;
    profile_rhob.open(input_rhob_filename.c_str());
    if (!profile_rhob.good()) {
        music_message << "Init::InitTJb: "
                      << "Can not open the initial file: "
                      << input_rhob_filename;
        music_message.flush("error");
        exit(1);
    }

    const int neta = arena_current.nEta();
    double temp_profile_ed[neta];
    double temp_profile_rhob[neta];

    double dummy;
    for (int ieta = 0; ieta < neta; ieta++) {
        profile_ed >> dummy >> temp_profile_ed[ieta];
        profile_rhob >> dummy >> temp_profile_rhob[ieta];
    }
    profile_ed.close();
    profile_rhob.close();

    const int nx = arena_current.nX();
    const int ny = arena_current.nY();
    for (int ieta = 0; ieta < neta; ieta++) {
        double rhob = temp_profile_rhob[ieta];
        double epsilon = temp_profile_ed[ieta]/hbarc;   // fm^-4
        for (int ix = 0; ix < nx; ix++) {
            for (int iy = 0; iy< ny; iy++) {
                // set all values in the grid element:
                arena_current(ix, iy, ieta).epsilon = epsilon;
                arena_current(ix, iy, ieta).rhob    = rhob;

                arena_current(ix, iy, ieta).u[0] = 1.0;
                arena_current(ix, iy, ieta).u[1] = 0.0;
                arena_current(ix, iy, ieta).u[2] = 0.0;
                arena_current(ix, iy, ieta).u[3] = 0.0;

                arena_prev(ix, iy, ieta) = arena_current(ix, iy, ieta);
            }
        }
    }
}

void Init::initial_IPGlasma_XY(int ieta, SCGrid &arena_prev,
                               SCGrid &arena_current) {
    ifstream profile(DATA.initName.c_str());

    std::string dummy;
    // read the information line
    std::getline(profile, dummy);

    const int nx = arena_current.nX();
    const int ny = arena_current.nY();
    double temp_profile_ed[nx][ny];
    double temp_profile_utau[nx][ny];
    double temp_profile_ux[nx][ny];
    double temp_profile_uy[nx][ny];

    // read the one slice
    double density, dummy1, dummy2, dummy3;
    double ux, uy, utau;
    for (int ix = 0; ix < nx; ix++) {
        for (int iy = 0; iy < ny; iy++) {
            profile >> dummy1 >> dummy2 >> dummy3
                    >> density >> utau >> ux >> uy
                    >> dummy  >> dummy  >> dummy  >> dummy;
            temp_profile_ed[ix][iy] = density;
            temp_profile_ux[ix][iy] = ux;
            temp_profile_uy[ix][iy] = uy;
            temp_profile_utau[ix][iy] = sqrt(1. + ux*ux + uy*uy);
            if (ix == 0 && iy == 0) {
                DATA.x_size = -dummy2*2;
                DATA.y_size = -dummy3*2;
                if (omp_get_thread_num() == 0) {
                    music_message << "eta_size=" << DATA.eta_size
                                  << ", x_size=" << DATA.x_size
                                  << ", y_size=" << DATA.y_size;
                    music_message.flush("info");
                }
            }
        }
    }
    profile.close();

    double eta = (DATA.delta_eta)*ieta - (DATA.eta_size)/2.0;
    double eta_envelop_ed = eta_profile_plateau(eta, DATA.eta_flat/2.0,
                                                DATA.eta_fall_off);
    int entropy_flag = DATA.initializeEntropy;
    for (int ix = 0; ix < nx; ix++) {
        for (int iy = 0; iy< ny; iy++) {
            double rhob = 0.0;
            double epsilon = 0.0;
            if (entropy_flag == 0) {
                epsilon = (temp_profile_ed[ix][iy]*eta_envelop_ed
                           *DATA.sFactor/hbarc);  // 1/fm^4
            } else {
                double local_sd = (temp_profile_ed[ix][iy]*DATA.sFactor
                                   *eta_envelop_ed);
                epsilon = eos.get_s2e(local_sd, rhob);
            }
            epsilon = std::max(Util::small_eps, epsilon);

            arena_current(ix, iy, ieta).epsilon = epsilon;
            arena_current(ix, iy, ieta).rhob = rhob;

            arena_current(ix, iy, ieta).u[0] = temp_profile_utau[ix][iy];
            arena_current(ix, iy, ieta).u[1] = temp_profile_ux[ix][iy];
            arena_current(ix, iy, ieta).u[2] = temp_profile_uy[ix][iy];
            arena_current(ix, iy, ieta).u[3] = 0.0;

            arena_prev(ix, iy, ieta) = arena_current(ix, iy, ieta);
        }
    }
}

void Init::initial_IPGlasma_XY_with_pi(int ieta, SCGrid &arena_prev,
                                       SCGrid &arena_current) {
    // Initial_profile == 9 : full T^\mu\nu
    // Initial_profile == 91: e and u^\mu
    // Initial_profile == 92: e only
    // Initial_profile == 93: e, u^\mu, and pi^\mu\nu, no bulk Pi
    double tau0 = DATA.tau0;
    ifstream profile(DATA.initName.c_str());

    std::string dummy;
    // read the information line
    std::getline(profile, dummy);

    const int nx = arena_current.nX();
    const int ny = arena_current.nY();
    std::vector<double> temp_profile_ed(nx*ny, 0.0);
    std::vector<double> temp_profile_utau(nx*ny, 0.0);
    std::vector<double> temp_profile_ux(nx*ny, 0.0);
    std::vector<double> temp_profile_uy(nx*ny, 0.0);
    std::vector<double> temp_profile_ueta(nx*ny, 0.0);
    std::vector<double> temp_profile_pitautau(nx*ny, 0.0);
    std::vector<double> temp_profile_pitaux(nx*ny, 0.0);
    std::vector<double> temp_profile_pitauy(nx*ny, 0.0);
    std::vector<double> temp_profile_pitaueta(nx*ny, 0.0);
    std::vector<double> temp_profile_pixx(nx*ny, 0.0);
    std::vector<double> temp_profile_pixy(nx*ny, 0.0);
    std::vector<double> temp_profile_pixeta(nx*ny, 0.0);
    std::vector<double> temp_profile_piyy(nx*ny, 0.0);
    std::vector<double> temp_profile_piyeta(nx*ny, 0.0);
    std::vector<double> temp_profile_pietaeta(nx*ny, 0.0);

    // read the one slice
    double density, dummy1, dummy2, dummy3;
    double ux, uy, utau, ueta;
    double pitautau, pitaux, pitauy, pitaueta;
    double pixx, pixy, pixeta, piyy, piyeta, pietaeta;
    for (int ix = 0; ix < nx; ix++) {
        for (int iy = 0; iy < ny; iy++) {
            int idx = iy + ix*ny;
            std::getline(profile, dummy);
            std::stringstream ss(dummy);
            ss >> dummy1 >> dummy2 >> dummy3
               >> density >> utau >> ux >> uy >> ueta
               >> pitautau >> pitaux >> pitauy >> pitaueta
               >> pixx >> pixy >> pixeta >> piyy >> piyeta >> pietaeta;
            ueta = ueta*tau0;
            temp_profile_ed    [idx] = density;
            temp_profile_ux    [idx] = ux;
            temp_profile_uy    [idx] = uy;
            temp_profile_ueta  [idx] = ueta;
            temp_profile_utau  [idx] = sqrt(1. + ux*ux + uy*uy + ueta*ueta);
            temp_profile_pixx  [idx] = pixx*DATA.sFactor;
            temp_profile_pixy  [idx] = pixy*DATA.sFactor;
            temp_profile_pixeta[idx] = pixeta*tau0*DATA.sFactor;
            temp_profile_piyy  [idx] = piyy*DATA.sFactor;
            temp_profile_piyeta[idx] = piyeta*tau0*DATA.sFactor;

            utau = temp_profile_utau[idx];
            temp_profile_pietaeta[idx] = (
                (2.*(  ux*uy*temp_profile_pixy[idx]
                     + ux*ueta*temp_profile_pixeta[idx]
                     + uy*ueta*temp_profile_piyeta[idx])
                 - (utau*utau - ux*ux)*temp_profile_pixx[idx]
                 - (utau*utau - uy*uy)*temp_profile_piyy[idx])
                /(utau*utau - ueta*ueta));
            temp_profile_pitaux  [idx] = (1./utau
                *(  temp_profile_pixx[idx]*ux
                  + temp_profile_pixy[idx]*uy
                  + temp_profile_pixeta[idx]*ueta));
            temp_profile_pitauy  [idx] = (1./utau
                *(  temp_profile_pixy[idx]*ux
                  + temp_profile_piyy[idx]*uy
                  + temp_profile_piyeta[idx]*ueta));
            temp_profile_pitaueta[idx] = (1./utau
                *(  temp_profile_pixeta[idx]*ux
                  + temp_profile_piyeta[idx]*uy
                  + temp_profile_pietaeta[idx]*ueta));
            temp_profile_pitautau[idx] = (1./utau
                *(  temp_profile_pitaux[idx]*ux
                  + temp_profile_pitauy[idx]*uy
                  + temp_profile_pitaueta[idx]*ueta));
            if (ix == 0 && iy == 0) {
                DATA.x_size = -dummy2*2;
                DATA.y_size = -dummy3*2;
                if (omp_get_thread_num() == 0) {
                    music_message << "eta_size=" << DATA.eta_size
                                  << ", x_size=" << DATA.x_size
                                  << ", y_size=" << DATA.y_size;
                    music_message.flush("info");
                }
            }
        }
    }
    profile.close();

    double eta = (DATA.delta_eta)*(ieta) - (DATA.eta_size)/2.0;
    double eta_envelop_ed = eta_profile_plateau(eta, DATA.eta_flat/2.0,
                                                DATA.eta_fall_off);
    int entropy_flag = DATA.initializeEntropy;
    for (int ix = 0; ix < nx; ix++) {
        for (int iy = 0; iy< ny; iy++) {
            int idx = iy + ix*ny;
            double rhob = 0.0;
            double epsilon = 0.0;
            if (entropy_flag == 0) {
                epsilon = (temp_profile_ed[idx]*eta_envelop_ed
                           *DATA.sFactor/hbarc);  // 1/fm^4
            } else {
                double local_sd = (temp_profile_ed[idx]*DATA.sFactor
                                   *eta_envelop_ed);
                epsilon = eos.get_s2e(local_sd, rhob);
            }
            epsilon = std::max(Util::small_eps, epsilon);

            arena_current(ix, iy, ieta).epsilon = epsilon;
            arena_current(ix, iy, ieta).rhob = rhob;

            if (DATA.Initial_profile == 92) {
                arena_current(ix, iy, ieta).u[0] = 1.0;
                arena_current(ix, iy, ieta).u[1] = 0.0;
                arena_current(ix, iy, ieta).u[2] = 0.0;
                arena_current(ix, iy, ieta).u[3] = 0.0;
            } else {
                arena_current(ix, iy, ieta).u[0] = temp_profile_utau[idx];
                arena_current(ix, iy, ieta).u[1] = temp_profile_ux[idx];
                arena_current(ix, iy, ieta).u[2] = temp_profile_uy[idx];
                arena_current(ix, iy, ieta).u[3] = temp_profile_ueta[idx];
            }

            if (DATA.Initial_profile == 9 || DATA.Initial_profile == 93) {
                arena_current(ix, iy, ieta).Wmunu[0] = temp_profile_pitautau[idx];
                arena_current(ix, iy, ieta).Wmunu[1] = temp_profile_pitaux[idx];
                arena_current(ix, iy, ieta).Wmunu[2] = temp_profile_pitauy[idx];
                arena_current(ix, iy, ieta).Wmunu[3] = temp_profile_pitaueta[idx];
                arena_current(ix, iy, ieta).Wmunu[4] = temp_profile_pixx[idx];
                arena_current(ix, iy, ieta).Wmunu[5] = temp_profile_pixy[idx];
                arena_current(ix, iy, ieta).Wmunu[6] = temp_profile_pixeta[idx];
                arena_current(ix, iy, ieta).Wmunu[7] = temp_profile_piyy[idx];
                arena_current(ix, iy, ieta).Wmunu[8] = temp_profile_piyeta[idx];
                arena_current(ix, iy, ieta).Wmunu[9] = temp_profile_pietaeta[idx];

                if (DATA.Initial_profile == 9) {
                    double pressure = eos.get_pressure(epsilon, rhob);
                    arena_current(ix, iy, ieta).pi_b = epsilon/3. - pressure;
                }
            }
            arena_prev(ix, iy, ieta) = arena_current(ix, iy, ieta);
        }
    }
}

void Init::initial_ROOT_XYEta_with_pi(SCGrid &arena_prev, SCGrid &arena_current) {
    #ifdef ROOT_FOUND
    // Initial_profile == 94 : full T^\mu\nu(x, y, eta)
    double tau0 = DATA.tau0;
    //ifstream profile(DATA.initName.c_str());

    const int nx = arena_current.nX();
    const int ny = arena_current.nY();
    const int neta = arena_current.nEta();

    //Open file and get histograms
    TFile* fIC = TFile::Open(DATA.initName.data(),"READ");
    TH3D* hE = (TH3D*) fIC->Get("energy_density");
    TH3D* hux = (TH3D*) fIC->Get("ux");
    TH3D* huy = (TH3D*) fIC->Get("uy");
    TH3D* hueta = (TH3D*) fIC->Get("ueta");
    TH3D* hpitautau = (TH3D*) fIC->Get("pi_tt");
    TH3D* hpitaux = (TH3D*) fIC->Get("pi_tx");
    TH3D* hpitauy = (TH3D*) fIC->Get("pi_ty");
    TH3D* hpitaueta = (TH3D*) fIC->Get("pi_teta");
    TH3D* hpixx = (TH3D*) fIC->Get("pi_xx");
    TH3D* hpixy = (TH3D*) fIC->Get("pi_xy");
    TH3D* hpixeta = (TH3D*) fIC->Get("pi_xeta");
    TH3D* hpiyy = (TH3D*) fIC->Get("pi_yy");
    TH3D* hpiyeta = (TH3D*) fIC->Get("pi_yeta");
    TH3D* hpietaeta = (TH3D*) fIC->Get("pi_etaeta");


    for (int ieta=0; ieta<neta;++ieta){

        std::vector<double> temp_profile_ed(nx*ny, 0.0);
        std::vector<double> temp_profile_utau(nx*ny, 0.0);
        std::vector<double> temp_profile_ux(nx*ny, 0.0);
        std::vector<double> temp_profile_uy(nx*ny, 0.0);
        std::vector<double> temp_profile_ueta(nx*ny, 0.0);
        std::vector<double> temp_profile_pitautau(nx*ny, 0.0);
        std::vector<double> temp_profile_pitaux(nx*ny, 0.0);
        std::vector<double> temp_profile_pitauy(nx*ny, 0.0);
        std::vector<double> temp_profile_pitaueta(nx*ny, 0.0);
        std::vector<double> temp_profile_pixx(nx*ny, 0.0);
        std::vector<double> temp_profile_pixy(nx*ny, 0.0);
        std::vector<double> temp_profile_pixeta(nx*ny, 0.0);
        std::vector<double> temp_profile_piyy(nx*ny, 0.0);
        std::vector<double> temp_profile_piyeta(nx*ny, 0.0);
        std::vector<double> temp_profile_pietaeta(nx*ny, 0.0);

        // read one slice
        double density;
        double ux, uy, utau, ueta;
        double pitautau, pitaux, pitauy, pitaueta;
        double pixx, pixy, pixeta, piyy, piyeta, pietaeta;
        for (int ix = 0; ix < nx; ix++) {
            for (int iy = 0; iy < ny; iy++) {
                int idx = iy + ix*ny;

                density = hE->GetBinContent(ix+1,iy+1,ieta+1);
                ux = hux->GetBinContent(ix+1,iy+1,ieta+1);
                uy = huy->GetBinContent(ix+1,iy+1,ieta+1);
                ueta = hueta->GetBinContent(ix+1,iy+1,ieta+1);
                pitautau = hpitautau->GetBinContent(ix+1,iy+1,ieta+1);
                pitaux = hpitaux->GetBinContent(ix+1,iy+1,ieta+1);
                pitauy = hpitauy->GetBinContent(ix+1,iy+1,ieta+1);
                pitaueta = hpitaueta->GetBinContent(ix+1,iy+1,ieta+1);
                pixx = hpixx->GetBinContent(ix+1,iy+1,ieta+1);
                pixy = hpixy->GetBinContent(ix+1,iy+1,ieta+1);
                pixeta = hpixeta->GetBinContent(ix+1,iy+1,ieta+1);
                piyy = hpiyy->GetBinContent(ix+1,iy+1,ieta+1);
                piyeta = hpiyeta->GetBinContent(ix+1,iy+1,ieta+1);
                pietaeta = hpietaeta->GetBinContent(ix+1,iy+1,ieta+1);

                ueta = ueta*tau0;
                temp_profile_ed    [idx] = density;
                temp_profile_ux    [idx] = ux;
                temp_profile_uy    [idx] = uy;
                temp_profile_ueta  [idx] = ueta;
                temp_profile_utau  [idx] = sqrt(1. + ux*ux + uy*uy + ueta*ueta);
                temp_profile_pixx  [idx] = pixx*DATA.sFactor;
                temp_profile_pixy  [idx] = pixy*DATA.sFactor;
                temp_profile_pixeta[idx] = pixeta*tau0*DATA.sFactor;
                temp_profile_piyy  [idx] = piyy*DATA.sFactor;
                temp_profile_piyeta[idx] = piyeta*tau0*DATA.sFactor;

                utau = temp_profile_utau[idx];
                temp_profile_pietaeta[idx] = (
                    (2.*(  ux*uy*temp_profile_pixy[idx]
                         + ux*ueta*temp_profile_pixeta[idx]
                         + uy*ueta*temp_profile_piyeta[idx])
                     - (utau*utau - ux*ux)*temp_profile_pixx[idx]
                     - (utau*utau - uy*uy)*temp_profile_piyy[idx])
                    /(utau*utau - ueta*ueta));
                temp_profile_pitaux  [idx] = (1./utau
                    *(  temp_profile_pixx[idx]*ux
                      + temp_profile_pixy[idx]*uy
                      + temp_profile_pixeta[idx]*ueta));
                temp_profile_pitauy  [idx] = (1./utau
                    *(  temp_profile_pixy[idx]*ux
                      + temp_profile_piyy[idx]*uy
                      + temp_profile_piyeta[idx]*ueta));
                temp_profile_pitaueta[idx] = (1./utau
                    *(  temp_profile_pixeta[idx]*ux
                      + temp_profile_piyeta[idx]*uy
                      + temp_profile_pietaeta[idx]*ueta));
                temp_profile_pitautau[idx] = (1./utau
                    *(  temp_profile_pitaux[idx]*ux
                      + temp_profile_pitauy[idx]*uy
                      + temp_profile_pitaueta[idx]*ueta));
            }
        }

        int entropy_flag = DATA.initializeEntropy;
        for (int ix = 0; ix < nx; ix++) {
            for (int iy = 0; iy< ny; iy++) {
                int idx = iy + ix*ny;
                double rhob = 0.0;
                double epsilon = 0.0;
                epsilon = (temp_profile_ed[idx]*DATA.sFactor/hbarc);  // 1/fm^4
                epsilon = std::max(Util::small_eps, epsilon);

                arena_current(ix, iy, ieta).epsilon = epsilon;
                arena_current(ix, iy, ieta).rhob = rhob;


                arena_current(ix, iy, ieta).u[0] = temp_profile_utau[idx];
                arena_current(ix, iy, ieta).u[1] = temp_profile_ux[idx];
                arena_current(ix, iy, ieta).u[2] = temp_profile_uy[idx];
                arena_current(ix, iy, ieta).u[3] = temp_profile_ueta[idx];

                arena_current(ix, iy, ieta).Wmunu[0] = temp_profile_pitautau[idx];
                arena_current(ix, iy, ieta).Wmunu[1] = temp_profile_pitaux[idx];
                arena_current(ix, iy, ieta).Wmunu[2] = temp_profile_pitauy[idx];
                arena_current(ix, iy, ieta).Wmunu[3] = temp_profile_pitaueta[idx];
                arena_current(ix, iy, ieta).Wmunu[4] = temp_profile_pixx[idx];
                arena_current(ix, iy, ieta).Wmunu[5] = temp_profile_pixy[idx];
                arena_current(ix, iy, ieta).Wmunu[6] = temp_profile_pixeta[idx];
                arena_current(ix, iy, ieta).Wmunu[7] = temp_profile_piyy[idx];
                arena_current(ix, iy, ieta).Wmunu[8] = temp_profile_piyeta[idx];
                arena_current(ix, iy, ieta).Wmunu[9] = temp_profile_pietaeta[idx];

                double pressure = eos.get_pressure(epsilon, rhob);
                if ( DATA.Initial_profile == 94) arena_current(ix, iy, ieta).pi_b = epsilon/3. - pressure;
                if ( DATA.Initial_profile == 95) arena_current(ix, iy, ieta).pi_b = -pressure;


                arena_prev(ix, iy, ieta) = arena_current(ix, iy, ieta);
            }
        }
    }
    fIC->Close();
    #else
    const int nx = arena_current.nX();
    const int ny = arena_current.nY();
    const int neta = arena_current.nEta();
    for (int ieta=0; ieta<neta;++ieta){
        for (int ix = 0; ix < nx; ix++) {
            for (int iy = 0; iy< ny; iy++) {
                arena_current(ix, iy, ieta).epsilon = 1E-12;
                arena_current(ix, iy, ieta).u[0] = 1.;
                arena_current(ix, iy, ieta).u[1] = .0;
                arena_current(ix, iy, ieta).u[2] = .0;
                arena_current(ix, iy, ieta).u[3] = .0;

                arena_current(ix, iy, ieta).Wmunu[0] = .0;
                arena_current(ix, iy, ieta).Wmunu[1] = .0;
                arena_current(ix, iy, ieta).Wmunu[2] = .0;
                arena_current(ix, iy, ieta).Wmunu[3] = .0;
                arena_current(ix, iy, ieta).Wmunu[4] = .0;
                arena_current(ix, iy, ieta).Wmunu[5] = .0;
                arena_current(ix, iy, ieta).Wmunu[6] = .0;
                arena_current(ix, iy, ieta).Wmunu[7] = .0;
                arena_current(ix, iy, ieta).Wmunu[8] = .0;
                arena_current(ix, iy, ieta).Wmunu[9] = .0;

                arena_current(ix, iy, ieta).pi_b = .0;
            }
        }
    }

    #endif
}

void Init::initial_MCGlb_with_rhob(SCGrid &arena_prev, SCGrid &arena_current) {
    // first load in the transverse profile
    ifstream profile_TA(DATA.initName_TA.c_str());
    ifstream profile_TB(DATA.initName_TB.c_str());

    const int nx = arena_current.nX();
    const int ny = arena_current.nY();
    const int neta = arena_current.nEta();
    double temp_profile_TA[nx][ny];
    double temp_profile_TB[nx][ny];
    double N_B = 0.0;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            profile_TA >> temp_profile_TA[i][j];
            profile_TB >> temp_profile_TB[i][j];
            N_B += temp_profile_TA[i][j] + temp_profile_TB[i][j];
        }
    }
    profile_TA.close();
    profile_TB.close();
    N_B *= DATA.delta_x*DATA.delta_y;
    double total_energy = DATA.ecm/2.*N_B;
    music_message << "sqrt{s} = " << DATA.ecm << " GeV, "
                  << "beam rapidity = " << DATA.beam_rapidity << ", "
                  << "total energy = " << total_energy << " GeV, "
                  << "N_B = " << N_B;
    music_message.flush("info");

    double T_tau_t = 0.0;
    #pragma omp parallel for reduction(+: T_tau_t)
    for (int ieta = 0; ieta < neta; ieta++) {
        double eta = (DATA.delta_eta)*ieta - (DATA.eta_size)/2.0;
        if (DATA.boost_invariant) {
            eta = 0.0;
        }
        double eta_rhob_left  = eta_rhob_left_factor(eta);
        double eta_rhob_right = eta_rhob_right_factor(eta);

        for (int ix = 0; ix < nx; ix++) {
            for (int iy = 0; iy< ny; iy++) {
                double rhob = 0.0;
                double epsilon = 0.0;
                if (DATA.turn_on_rhob == 1) {
                    rhob = (  temp_profile_TA[ix][iy]*eta_rhob_right
                            + temp_profile_TB[ix][iy]*eta_rhob_left);
                } else {
                    rhob = 0.0;
                }

                if (DATA.Initial_profile == 11) {
                    const double eta_0 = DATA.eta_flat/2.;
                    const double sigma_eta = DATA.eta_fall_off;
                    const double E_norm = energy_eta_profile_normalisation(
                                                0.0, eta_0, sigma_eta);
                    const double Pz_norm = Pz_eta_profile_normalisation(
                                                eta_0, sigma_eta);
                    const double norm_even = (
                            1./(DATA.tau0*E_norm)
                            *Util::m_N*cosh(DATA.beam_rapidity));
                    const double norm_odd = (
                            DATA.beam_rapidity/(DATA.tau0*Pz_norm)
                            *Util::m_N*sinh(DATA.beam_rapidity));
                    double eta_envelop = eta_profile_plateau(
                                                    eta, eta_0, sigma_eta);
                    epsilon = (
                        ((  (temp_profile_TA[ix][iy] + temp_profile_TB[ix][iy])
                           *norm_even
                          + (temp_profile_TA[ix][iy] - temp_profile_TB[ix][iy])
                            *norm_odd*eta/DATA.beam_rapidity)*eta_envelop)
                        /Util::hbarc);
                } else if (DATA.Initial_profile == 111) {
                    double y_CM = atanh(
                        (temp_profile_TA[ix][iy] - temp_profile_TB[ix][iy])
                        /(temp_profile_TA[ix][iy] + temp_profile_TB[ix][iy]
                          + Util::small_eps)
                        *tanh(DATA.beam_rapidity));
                    // local energy density [1/fm]
                    double E_lrf = (
                        (temp_profile_TA[ix][iy] + temp_profile_TB[ix][iy])
                        *Util::m_N*cosh(DATA.beam_rapidity)/Util::hbarc);
                    double eta0 = std::min(DATA.eta_flat/2.0,
                                    std::abs(DATA.beam_rapidity - y_CM));
                    double eta_envelop = eta_profile_plateau(
                                    eta - y_CM, eta0, DATA.eta_fall_off);
                    double E_norm = (
                        DATA.tau0*energy_eta_profile_normalisation(
                                    y_CM, eta0, DATA.eta_fall_off));
                    epsilon = E_lrf*eta_envelop/E_norm;
                }
                epsilon = std::max(Util::small_eps, epsilon);

                arena_current(ix, iy, ieta).epsilon = epsilon;
                arena_current(ix, iy, ieta).rhob = rhob;

                arena_current(ix, iy, ieta).u[0] = 1.0;
                arena_current(ix, iy, ieta).u[1] = 0.0;
                arena_current(ix, iy, ieta).u[2] = 0.0;
                arena_current(ix, iy, ieta).u[3] = 0.0;

                T_tau_t += epsilon*cosh(eta);
            }
        }
    }
    T_tau_t *= DATA.tau0*DATA.delta_eta*DATA.delta_x*DATA.delta_y*Util::hbarc;
    double norm = total_energy/T_tau_t;
    music_message << "energy norm = " << norm;
    music_message.flush("info");

    // renormalize the system's energy density
    #pragma omp parallel for collapse(3)
    for (int ieta = 0; ieta < neta; ieta++) {
        for (int ix = 0; ix < nx; ix++) {
            for (int iy = 0; iy< ny; iy++) {
                //arena_current(ix, iy, ieta).epsilon *= norm;
                arena_prev(ix, iy, ieta) = arena_current(ix, iy, ieta);
            }
        }
    }
}


void Init::initial_with_zero_XY(int ieta, SCGrid &arena_prev,
                                SCGrid &arena_current) {
    const int nx = arena_current.nX();
    const int ny = arena_current.nY();
    double u[4] = {1.0, 0.0, 0.0, 0.0};
    for (int ix = 0; ix < nx; ix++) {
        for (int iy = 0; iy < ny; iy++) {
            double rhob = 0.0;
            double epsilon = 1e-12;

            arena_current(ix, iy, ieta).epsilon = epsilon;
            arena_current(ix, iy, ieta).rhob = rhob;

            arena_current(ix, iy, ieta).u[0] = u[0];
            arena_current(ix, iy, ieta).u[1] = u[1];
            arena_current(ix, iy, ieta).u[2] = u[2];
            arena_current(ix, iy, ieta).u[3] = u[3];

            arena_prev(ix, iy, ieta) = arena_current(ix, iy, ieta);
        }
    }
}


void Init::initial_UMN_with_rhob(SCGrid &arena_prev, SCGrid &arena_current) {
    // first load in the transverse profile
    ifstream profile(DATA.initName.c_str());

    if (!profile) {
        music_message << "Can not open file: " << DATA.initName;
        music_message.flush("error");
        exit(1);
    }
    std::string dummy_s;
    std::getline(profile, dummy_s);

    const int nx   = arena_current.nX();
    const int ny   = arena_current.nY();
    const int neta = arena_current.nEta();

    double dummy;
    double ed_local, rhob_local;
    for (int ieta = 0; ieta < neta; ieta++) {
        for (int ix = 0; ix < nx; ix++) {
            for (int iy = 0; iy< ny; iy++) {
                profile >> dummy >> dummy >> dummy >> rhob_local >> ed_local;
                double rhob    = rhob_local;
                double epsilon = ed_local*DATA.sFactor/hbarc;    // 1/fm^4

                epsilon = std::max(Util::small_eps, epsilon);

                arena_current(ix, iy, ieta).epsilon = epsilon;
                arena_current(ix, iy, ieta).rhob = rhob;

                arena_current(ix, iy, ieta).u[0] = 1.0;
                arena_current(ix, iy, ieta).u[1] = 0.0;
                arena_current(ix, iy, ieta).u[2] = 0.0;
                arena_current(ix, iy, ieta).u[3] = 0.0;

                arena_prev(ix, iy, ieta) = arena_current(ix, iy, ieta);
            }
        }
    }
    profile.close();
}

void Init::initial_AMPT_XY(int ieta, SCGrid &arena_prev,
                           SCGrid &arena_current) {
    double u[4] = {1.0, 0.0, 0.0, 0.0};
    EnergyFlowVec j_mu = {0.0, 0.0, 0.0, 0.0};

    double eta = (DATA.delta_eta)*ieta - (DATA.eta_size)/2.0;
    double tau0 = DATA.tau0;
    const int nx = arena_current.nX();
    const int ny = arena_current.nY();
    for (int ix = 0; ix < nx; ix++) {
        double x_local = - DATA.x_size/2. + ix*DATA.delta_x;
        for (int iy = 0; iy < ny; iy++) {
            double y_local = - DATA.y_size/2. + iy*DATA.delta_y;
            double rhob = 0.0;
            double epsilon = 0.0;
            if (DATA.turn_on_rhob == 1) {
                rhob = hydro_source_terms_ptr.lock()->get_hydro_rhob_source_before_tau(
                                                tau0, x_local, y_local, eta);
            } else {
                rhob = 0.0;
            }

            hydro_source_terms_ptr.lock()->get_hydro_energy_source_before_tau(
                                    tau0, x_local, y_local, eta, j_mu);

            epsilon = j_mu[0];           // 1/fm^4

            epsilon = std::max(Util::small_eps, epsilon);

            arena_current(ix, iy, ieta).epsilon = epsilon;
            arena_current(ix, iy, ieta).rhob = rhob;

            arena_current(ix, iy, ieta).u[0] = u[0];
            arena_current(ix, iy, ieta).u[1] = u[1];
            arena_current(ix, iy, ieta).u[2] = u[2];
            arena_current(ix, iy, ieta).u[3] = u[3];

            arena_prev(ix, iy, ieta) = arena_current(ix, iy, ieta);
        }
    }
}


void Init::get_jetscape_preequilibrium_vectors(
        vector<double> e_in,
        vector<double> u_tau_in, vector<double> u_x_in,
        vector<double> u_y_in,   vector<double> u_eta_in,
        vector<double> pi_00_in, vector<double> pi_01_in,
        vector<double> pi_02_in, vector<double> pi_03_in,
        vector<double> pi_11_in, vector<double> pi_12_in,
        vector<double> pi_13_in, vector<double> pi_22_in,
        vector<double> pi_23_in, vector<double> pi_33_in,
        vector<double> Bulk_pi_in) {
    jetscape_initial_energy_density = e_in;
    jetscape_initial_u_tau          = u_tau_in;
    jetscape_initial_u_x            = u_x_in;
    jetscape_initial_u_y            = u_y_in;
    jetscape_initial_u_eta          = u_eta_in;
    jetscape_initial_pi_00          = pi_00_in;
    jetscape_initial_pi_01          = pi_01_in;
    jetscape_initial_pi_02          = pi_02_in;
    jetscape_initial_pi_03          = pi_03_in;
    jetscape_initial_pi_11          = pi_11_in;
    jetscape_initial_pi_12          = pi_12_in;
    jetscape_initial_pi_13          = pi_13_in;
    jetscape_initial_pi_22          = pi_22_in;
    jetscape_initial_pi_23          = pi_23_in;
    jetscape_initial_pi_33          = pi_33_in;
    jetscape_initial_bulk_pi        = Bulk_pi_in;
}


void Init::initial_with_jetscape(int ieta, SCGrid &arena_prev,
                                 SCGrid &arena_current) {
    const int nx = arena_current.nX();
    const int ny = arena_current.nY();
    //const int neta = arena_current.nEta();

    for (int ix = 0; ix < nx; ix++) {
        for (int iy = 0; iy< ny; iy++) {
            const double rhob = 0.0;
            double epsilon = 0.0;
            const int idx = iy + ix*ny + ieta*ny*nx;
            epsilon = (jetscape_initial_energy_density[idx]
                       *DATA.sFactor/hbarc);  // 1/fm^4
            epsilon = std::max(Util::small_eps, epsilon);

            arena_current(ix, iy, ieta).epsilon = epsilon;
            arena_current(ix, iy, ieta).rhob = rhob;

            arena_current(ix, iy, ieta).u[0] = jetscape_initial_u_tau[idx];
            arena_current(ix, iy, ieta).u[1] = jetscape_initial_u_x[idx];
            arena_current(ix, iy, ieta).u[2] = jetscape_initial_u_y[idx];
            arena_current(ix, iy, ieta).u[3] = DATA.tau0*jetscape_initial_u_eta[idx];

            arena_current(ix, iy, ieta).pi_b = jetscape_initial_bulk_pi[idx]/hbarc;

            arena_current(ix, iy, ieta).Wmunu[0] = jetscape_initial_pi_00[idx]/hbarc;
            arena_current(ix, iy, ieta).Wmunu[1] = jetscape_initial_pi_01[idx]/hbarc;
            arena_current(ix, iy, ieta).Wmunu[2] = jetscape_initial_pi_02[idx]/hbarc;
            arena_current(ix, iy, ieta).Wmunu[3] = jetscape_initial_pi_03[idx]/hbarc*DATA.tau0;
            arena_current(ix, iy, ieta).Wmunu[4] = jetscape_initial_pi_11[idx]/hbarc;
            arena_current(ix, iy, ieta).Wmunu[5] = jetscape_initial_pi_12[idx]/hbarc;
            arena_current(ix, iy, ieta).Wmunu[6] = jetscape_initial_pi_13[idx]/hbarc*DATA.tau0;
            arena_current(ix, iy, ieta).Wmunu[7] = jetscape_initial_pi_22[idx]/hbarc;
            arena_current(ix, iy, ieta).Wmunu[8] = jetscape_initial_pi_23[idx]/hbarc*DATA.tau0;
            arena_current(ix, iy, ieta).Wmunu[9] = jetscape_initial_pi_33[idx]/hbarc*DATA.tau0*DATA.tau0;

            arena_prev(ix, iy, ieta) = arena_current(ix, iy, ieta);
        }
    }
}

void Init::clean_up_jetscape_arrays() {
    // clean up
    jetscape_initial_energy_density.clear();
    jetscape_initial_u_tau.clear();
    jetscape_initial_u_x.clear();
    jetscape_initial_u_y.clear();
    jetscape_initial_u_eta.clear();
    jetscape_initial_pi_00.clear();
    jetscape_initial_pi_01.clear();
    jetscape_initial_pi_02.clear();
    jetscape_initial_pi_03.clear();
    jetscape_initial_pi_11.clear();
    jetscape_initial_pi_12.clear();
    jetscape_initial_pi_13.clear();
    jetscape_initial_pi_22.clear();
    jetscape_initial_pi_23.clear();
    jetscape_initial_pi_33.clear();
    jetscape_initial_bulk_pi.clear();
}


double Init::eta_profile_plateau(const double eta, const double eta_0,
                                 const double sigma_eta) const {
    // this function return the eta envelope profile for energy density
    // Hirano's plateau + Gaussian fall-off
    double res;
    double exparg1 = (std::abs(eta) - eta_0)/sigma_eta;
    double exparg = exparg1*exparg1/2.0;
    res = exp(-exparg*Util::theta(exparg1));
    return res;
}


double Init::energy_eta_profile_normalisation(
        const double y_CM, const double eta_0, const double sigma_eta) const {
    // this function returns the normalization of the eta envelope profile
    // for energy density
    //  \int deta eta_profile_plateau(eta - y_CM, eta_0, sigma_eta)*cosh(eta)
    double f1 = (  exp( eta_0)*erfc(-sqrt(0.5)*sigma_eta)
                 + exp(-eta_0)*erfc( sqrt(0.5)*sigma_eta));
    double f2 = sqrt(M_PI/2.)*sigma_eta*exp(sigma_eta*sigma_eta/2.);
    double f3 = sinh(eta_0 + y_CM) - sinh(-eta_0 + y_CM);
    double norm = cosh(y_CM)*f2*f1 + f3;
    return(norm);
}


double Init::Pz_eta_profile_normalisation(
        const double eta_0, const double sigma_eta) const {
    // this function returns the normalization of the eta envelope profile
    // for longitudinal momentum
    //  \int deta eta_profile_plateau(eta, eta_0, sigma_eta)*eta*sinh(eta)
    const double sigma_sq = sigma_eta*sigma_eta;
    double f1 = (  exp( eta_0)*(eta_0 + sigma_sq)*erfc(-sqrt(0.5)*sigma_eta)
                 - exp(-eta_0)*(eta_0 - sigma_sq)*erfc( sqrt(0.5)*sigma_eta));
    double f2 = sqrt(M_PI/2.)*sigma_eta*exp(sigma_sq/2.)/2.;
    double f3 = sigma_sq*sinh(eta_0);
    double f4 = 2.*eta_0*cosh(eta_0) - 2.*sinh(eta_0);
    double norm = 2.*(f2*f1 + f3) + f4;
    return(norm);
}

double Init::eta_profile_left_factor(const double eta) const {
    // this function return the eta envelope for projectile
    double res = eta_profile_plateau(
                    eta, DATA.eta_flat/2.0, DATA.eta_fall_off);
    if (std::abs(eta) < DATA.beam_rapidity) {
        res = (1. - eta/DATA.beam_rapidity)*res;
    } else {
        res = 0.0;
    }
    return(res);
}


double Init::eta_profile_right_factor(const double eta) const {
    // this function return the eta envelope for target
    double res = eta_profile_plateau(
                    eta, DATA.eta_flat/2.0, DATA.eta_fall_off);
    if (std::abs(eta) < DATA.beam_rapidity) {
        res = (1. + eta/DATA.beam_rapidity)*res;
    } else {
        res = 0.0;
    }
    return(res);
}

double Init::eta_rhob_profile_normalisation(const double eta) const {
    // this function return the eta envelope profile for net baryon density
    double res = 0.0;
    int profile_flag = DATA.initial_eta_rhob_profile;
    double eta_0 = DATA.eta_rhob_0;
    double tau0 = DATA.tau0;
    if (profile_flag == 1) {
        const double eta_width = DATA.eta_rhob_width;
        const double norm      = 1./(2.*sqrt(2*M_PI)*eta_width*tau0);
        const double exparg1   = (eta - eta_0)/eta_width;
        const double exparg2   = (eta + eta_0)/eta_width;
        res = norm*(exp(-exparg1*exparg1/2.0) + exp(-exparg2*exparg2/2.0));
    } else if (profile_flag == 2) {
        double eta_abs     = fabs(eta);
        double delta_eta_1 = DATA.eta_rhob_width_1;
        double delta_eta_2 = DATA.eta_rhob_width_2;
        double A           = DATA.eta_rhob_plateau_height;
        double exparg1     = (eta_abs - eta_0)/delta_eta_1;
        double exparg2     = (eta_abs - eta_0)/delta_eta_2;
        double theta;
        double norm = 1./(tau0*(sqrt(2.*M_PI)*delta_eta_1
                          + (1. - A)*sqrt(2.*M_PI)*delta_eta_2 + 2.*A*eta_0));
        if (eta_abs > eta_0)
            theta = 1.0;
        else
            theta = 0.0;
        res = norm*(theta*exp(-exparg1*exparg1/2.)
                    + (1. - theta)*(A + (1. - A)*exp(-exparg2*exparg2/2.)));
    }
    return res;
}


double Init::eta_rhob_left_factor(const double eta) const {
    double eta_0       = -std::abs(DATA.eta_rhob_0);
    double tau0        = DATA.tau0;
    double delta_eta_1 = DATA.eta_rhob_width_1;
    double delta_eta_2 = DATA.eta_rhob_width_2;
    double norm        = 2./(sqrt(M_PI)*tau0*(delta_eta_1 + delta_eta_2));
    double exp_arg     = 0.0;
    if (eta < eta_0) {
        exp_arg = (eta - eta_0)/delta_eta_1;
    } else {
        exp_arg = (eta - eta_0)/delta_eta_2;
    }
    double res = norm*exp(-exp_arg*exp_arg);
    return(res);
}


double Init::eta_rhob_right_factor(const double eta) const {
    double eta_0       = std::abs(DATA.eta_rhob_0);
    double tau0        = DATA.tau0;
    double delta_eta_1 = DATA.eta_rhob_width_1;
    double delta_eta_2 = DATA.eta_rhob_width_2;
    double norm        = 2./(sqrt(M_PI)*tau0*(delta_eta_1 + delta_eta_2));
    double exp_arg     = 0.0;
    if (eta < eta_0) {
        exp_arg = (eta - eta_0)/delta_eta_2;
    } else {
        exp_arg = (eta - eta_0)/delta_eta_1;
    }
    double res = norm*exp(-exp_arg*exp_arg);
    return(res);
}

void Init::output_initial_density_profiles(SCGrid &arena) {
    // this function outputs the 3d initial energy density profile
    // and net baryon density profile (if turn_on_rhob == 1)
    // for checking purpose
    music_message.info("output initial density profiles into a file... ");
    std::ofstream of("check_initial_density_profiles.dat");
    of << "# x(fm)  y(fm)  eta  ed(GeV/fm^3)";
    if (DATA.turn_on_rhob == 1)
        of << "  rhob(1/fm^3)";
    of << std::endl;
    for (int ieta = 0; ieta < arena.nEta(); ieta++) {
        double eta_local = (DATA.delta_eta)*ieta - (DATA.eta_size)/2.0;
        for(int ix = 0; ix < arena.nX(); ix++) {
            double x_local = -DATA.x_size/2. + ix*DATA.delta_x;
            for(int iy = 0; iy < arena.nY(); iy++) {
                double y_local = -DATA.y_size/2. + iy*DATA.delta_y;
                of << std::scientific << std::setw(18) << std::setprecision(8)
                   << x_local << "   " << y_local << "   "
                   << eta_local << "   " << arena(ix,iy,ieta).epsilon*hbarc;
                if (DATA.turn_on_rhob == 1) {
                    of << "   " << arena(ix,iy,ieta).rhob;
                }
                of << std::endl;
            }
        }
    }
    music_message.info("done!");
}

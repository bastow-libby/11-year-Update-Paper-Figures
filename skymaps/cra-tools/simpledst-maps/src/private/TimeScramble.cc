#include <SimpleDST.h>
#include <config.h>

#include <TChain.h>
#include <TH1D.h>
#include <TMath.h>
#include <TRandom.h>
#include <TStopwatch.h>

#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <memory>

#include <healpix_cxx/fitshandle.h>
#include <healpix_cxx/healpix_map.h>
#include <healpix_cxx/healpix_map_fitsio.h>
#include <healpix_cxx/pointing.h>

#include <photospline/splinetable.h>
#include <photospline/bspline.h>

#include <boost/program_options.hpp>

#include <astro/astro.h>
#include <astro/time.h>
#include <Direction.h>
#include <solardipole.h>

using namespace std;
namespace po = boost::program_options;

const double hour = 1 / 24.;
const double second = hour / 3600.;


typedef Healpix_Map<float> SkyMap;
typedef boost::shared_ptr<SkyMap> SkyMapPtr;


// Default values for inputs
struct ProgramOptions {

  // Input/output
  std::vector< std::string > input;
  std::string outdir = "./sample";
  std::string outfile = "CR_TS";

  // Time-scrambling inputs
  std::string detector = "IC";
  std::string yyyymmdd;
  int nInt;
  std::string method = "sid";
  int resample = 20;

  // Energy reconstruction
  std::string spline;
  std::vector<float> ebins;
  std::vector<float> sbins;

  // Solar dipole flags
  bool sd = false;
  bool sd2 = false;

};


// NOTES:
// - energy cuts should be rewritten as map placement
// - figure out warnings and throw them when...
//    - nInt is set but not yyyymmdd (could lead to weirdly-sized windows)
//    - ...
// - something's wrong with teh time scrambling checks. Running on the root
// file for IceTop gave me a 1 hour integration time and multiple time
// scrambling windows (should just be one without a specified nInt...)

void tScramble(po::variables_map vm, const ProgramOptions& opts);
int ITs125Cut(po::variables_map vm, SimpleDST dst, std::vector<float> sbins);
int ICenergyCut(po::variables_map vm, SimpleDST dst, photospline::splinetable<> &table, double zenith, std::vector<float> ebins);
double time_standard_ra(double mjd, Equatorial eq, std::string time_standard);


int main(int argc, char* argv[]) {

  ProgramOptions opts;

  po::options_description desc("Allowed options");
  desc.add_options()
      // Options used for all configurations
      ("help", "Produce help message")
      ("input", po::value<std::vector< std::string > >(&opts.input)
          ->multitoken()
          ->required(),
          "Input files")
      ("outdir", po::value<std::string>(&opts.outdir)
          ->default_value(opts.outdir),
          "Directory of output")
      ("outfile", po::value<std::string>(&opts.outfile)
          ->default_value(opts.outfile),
          "Base name for output file")
      ("yyyymmdd", po::value<std::string>(&opts.yyyymmdd),
          "Restrict scrambling to times within target date")
      ("detector", po::value<std::string>(&opts.detector)
          ->default_value(opts.detector),
          "Detector configuration (IC or IT)")
      ("method", po::value<std::string>(&opts.method)
          ->default_value(opts.method),
          "Time standard (sid|anti|solar|ext)")

      // Less-common options
      ("nInt", po::value<int>(&opts.nInt),
          "Integration time in hours")
      ("resample", po::value<int>(&opts.resample)
          ->default_value(opts.resample),
          "Number of times to resample background values (default: 20)")

      // Solar dipole flags
      ("sd", po::value<bool>(&opts.sd)
          ->default_value(opts.sd),
          "Correct for solar dipole for each event")
      ("sd2", po::value<bool>(&opts.sd2)
          ->default_value(opts.sd2),
          "Apply 2nd-order dipole correction")

      // IceCube specific options
      ("spline", po::value<std::string>(&opts.spline),
          "File containing spline tables")
      ("ebins", po::value< std::vector<float> >(&opts.ebins)
          ->multitoken(),
          "Energy bins")

      // Icetop specific options
      ("sbins", po::value< std::vector<float> >(&opts.sbins)
          ->multitoken(),
          "S125 bins")
  ;

  po::variables_map vm;

  // Disable short flags ('-') to allow use of negative signs in input
  po::store(po::parse_command_line(argc, argv, desc,
      po::command_line_style::unix_style ^ po::command_line_style::allow_short),
      vm);

  po::notify(vm);

  if (vm.count("help")) {
    std::cout << desc << "\n";
    return 1;
  }

  std::cout << "Input files are: " << endl;
  for (unsigned i = 0; i < opts.input.size(); i++)
      std::cout << opts.input[i] << endl;

  if (vm.count("ebins")) {
    std::cout << "Ebin values:" << endl;
    for (float val : opts.ebins)
        std::cout << val << " ";
    std::cout << endl;
  }

  if (vm.count("sbins")) {
    std::cout << "Sbin values:" << endl;
    for (float val : opts.sbins)
        std::cout << val << " ";
    std::cout << endl;
  }

  tScramble(vm, opts);

  return 0;
}


void tScramble(po::variables_map vm, const ProgramOptions& opts) {

  std::cout << "good start..." << endl;

  TStopwatch timer;
  timer.Start();

  // Read in spline tables if provided
  photospline::splinetable<> spline;
  if (vm.count("spline")) {
    std::cout << opts.spline << endl;
    spline.read_fits(opts.spline.c_str());
  }

  // Allow for a map for each energy (or S125) bin
  unsigned nMaps = 1;
  if (vm.count("ebins"))
    nMaps = opts.ebins.size() - 1;
  if (vm.count("sbins"))
    nMaps = opts.sbins.size() - 1;
  std::cout << "Number of maps: " << nMaps << endl;


  // Create map storage using a shared pointer
  const int NSide = 64;

  std::vector<SkyMapPtr> LocalMap;
  std::vector<SkyMapPtr> DataMap;
  std::vector<SkyMapPtr> BGMap;

  for (unsigned int i=0; i<nMaps; i++) {
    SkyMapPtr ebinMap(new SkyMap);
    ebinMap->SetNside(NSide, RING);
    ebinMap->fill(0.);
    LocalMap.push_back(ebinMap);
    DataMap.push_back(ebinMap);
    BGMap.push_back(ebinMap);
  }

  const char* masterTree;
  if (opts.detector == "IC")
    masterTree = "CutDST";
  else if (opts.detector == "IT")
    masterTree = "MasterTree";

  // Initialize the chain and read data
  //TChain* cutDST = new TChain(masterTree);
  auto cutDST = std::make_unique<TChain>(masterTree);
  for (unsigned i = 0; i < opts.input.size(); i++) {
    cutDST->Add(opts.input[i].c_str());
  }
  SimpleDST dst(cutDST.get(), opts.detector);
  std::cout << "Number of chained files: " << cutDST->GetNtrees() << endl;

  Long64_t nEntries = cutDST->GetEntries();
  std::cout << "Number of entries: " << nEntries << endl;

  // NOTE: needs to fail if nInt specified and not yyyymmdd
  // (could lead to weird short-window rounding errors)
  if ( vm.count("nInt") && !vm.count("yyyymmdd") ) {
    std::cout << "Specifying nInt requires a yyyymmdd argument" << endl;
    std::cout << "Raise a system error and exit" << endl;
  }

  // Use fixed start and stop times if yyyymmdd option provided
  double start_mjd, stop_mjd;
  if (vm.count("yyyymmdd")) {
    int yy = atoi(opts.yyyymmdd.substr(0, 4).c_str());
    int mm = atoi(opts.yyyymmdd.substr(5, 2).c_str());
    int dd = atoi(opts.yyyymmdd.substr(8, 2).c_str());
    astro::Time t(yy, mm, dd, 0, 0, 0);
    start_mjd = t.GetMJD();
    stop_mjd = start_mjd + 1;
  }

  // Otherwise, set start and stop time based on data
  else {
    cutDST->GetEntry(0);
    start_mjd = dst.ModJulDay;
    cutDST->GetEntry(nEntries - 1);
    stop_mjd = dst.ModJulDay;
  }

  // MJD1 will update to a new start time if scrambling window is passed
  double mjd1 = start_mjd;

  // Establish time-integration window
  int dt_hrs = opts.nInt;
  if (!vm.count("nInt")) {
    if (vm.count("yyyymmdd"))
      dt_hrs = 24;
    else
      dt_hrs = static_cast<int>(std::ceil(stop_mjd - start_mjd));
  }

  std::vector<Long64_t> nEvents(nMaps, 0);
  std::vector<Long64_t> nUsedEvents(nMaps, 0);

  const double alpha = 1. / opts.resample;

  // Integration time
  const double dt = dt_hrs * hour;
  std::cout << "Integration time = " << dt_hrs << " (hours) " << dt << " (day)\n";
  std::cout << "Reading " << nEntries << " entries...\n";

  // Setup histograms for storing time information
  std::vector<TH1D*> histMJD(nMaps);
  for (unsigned i = 0; i < nMaps; i++) {
    std::string histName = "histMJD_" + std::to_string(i);
    histMJD[i] = new TH1D(histName.c_str(), ";modified julian day;events",
        Int_t((stop_mjd - start_mjd) / (10. * second)), start_mjd, stop_mjd);
  }

  // Track the local coordinates
  std::vector< std::vector<Float_t> > LocCoord_theta(nMaps);
  std::vector< std::vector<Float_t> > LocCoord_phi(nMaps);

  // Timers to figure out what's taking so long...
  TStopwatch timer1;
  timer1.Start();

  //=====================================================================//
  // Begin iterating through events
  //=====================================================================//
  int dayCounter = 0;
  int validCounter = 0;
  int nevent = 0;
  double mjd = 0;

  // Read all events (potentially up to stop time)
  while ( nevent < nEntries && mjd < stop_mjd ) {

    cutDST->GetEntry(nevent);

    //==================================================================//
    // Basic tracking output
    //==================================================================//

    // Events before start of time window
    if (dst.ModJulDay < start_mjd) {
      if (nevent % 10000000 == 0)
        std::cout << "Processed " << nevent << " entries before starting..." << endl;
      nevent++;
      continue;
    }

    // Time to hit first entry
    dayCounter += 1;
    if (dayCounter == 1) {
      std::cout << "First entry: " << nevent << endl;
      timer1.Stop();
      printf("Time to first entry: %7.3fs\n", timer.RealTime());
    }

    // Counted event tracker
    if (dayCounter % 1000000 == 0) {
      std::cout << "Processed " << dayCounter << " entries in the right day of " << nevent+1 << " total entries..." << endl;
    }

    //==================================================================//
    // Event cuts
    //==================================================================//

    bool event_passed = true;

    // Extract zenith, azimuth, and fit_status information
    double zenith, azimuth;
    bool fitPassed;
    if (opts.detector == "IC") {
      zenith = dst.LLHZenithDeg * M_PI/180.;
      azimuth = dst.LLHAzimuthDeg * M_PI/180.;
      fitPassed = dst.isReco;
    }
    else if (opts.detector == "IT") {
      zenith = dst.SPZenith;
      azimuth = dst.SPAzimuth;
      fitPassed = (dst.SPFitStatus == 0);
    }

    // Throw away reconstructions too close to poles
    const float zLo = 0.002;           // 0.11 degrees
    const float zHi = M_PI - 0.002;    // 179.89 degrees
    if (zenith < zLo || zenith > zHi)
      event_passed = false;

    // Reconstruction cuts
    if (!fitPassed || std::isnan(zenith) || std::isnan(azimuth))
      event_passed = false;

    // Energy cuts for IceTop and IceCube
    int mapIdx = 0;
    if (vm.count("spline"))
      mapIdx = ICenergyCut(vm, dst, spline, zenith, opts.ebins);
    if (vm.count("sbins"))
      mapIdx = ITs125Cut(vm, dst, opts.sbins);
    // Energy bin of -1 is outside range
    if (mapIdx == -1)
      event_passed = false;


    mjd = dst.ModJulDay;
    if ( event_passed && mjd <= (mjd1+dt) ) {

      validCounter += 1;

      // Store local coordinates
      ++nEvents[mapIdx];
      LocCoord_theta[mapIdx].push_back(zenith);
      LocCoord_phi[mapIdx].push_back(azimuth);

      // Calculate coordinates in equatoral time
      Direction dir(zenith,azimuth);
      Equatorial eq = GetEquatorialFromDirection(dir, mjd);

      // RA can change depending on time standard (anti, ext, solar, sid)
      double ra = time_standard_ra(mjd, eq, opts.method);

      // Calculate solar dipole weighting
      double eventweight = 1.0;
      if (opts.sd)
        eventweight = solar_dipole(mjd, eq.ra, eq.dec, opts.sd2);
      pointing localDir(zenith, azimuth);
      SkyMapPtr tmp_map = LocalMap[mapIdx];
      int pixelID = tmp_map->ang2pix(localDir);
      (*tmp_map)[pixelID] += eventweight;

      // Write to map
      pointing eqDir(M_PI/2.-eq.dec, ra);
      tmp_map = DataMap[mapIdx];
      pixelID = tmp_map->ang2pix(eqDir);
      (*tmp_map)[pixelID] += eventweight;

      // Store time
      histMJD[mapIdx]->Fill(dst.ModJulDay);
    }

    // Need to check for time-scrambling behavior here
    if (mjd > (mjd1+dt) || mjd > stop_mjd || nevent+1 == nEntries) {

      for (unsigned mEntry = 0; mEntry<nMaps; mEntry++) {

        if (vm.count("ebins")) {
          std::cout << "Working on energy bin " << opts.ebins[mEntry] << "-"
               << opts.ebins[mEntry+1] << "GeV..." << endl;
        }
        if (vm.count("sbins")) {
          std::cout << "Working on s125 bin " << opts.sbins[mEntry] << " to "
               << opts.sbins[mEntry+1] << "s125..." << endl;
        }
        nUsedEvents[mEntry] += (nEvents[mEntry]);

        // Scramble the time
        std::cout << "  Scrambling time for (" << opts.resample << " x "
             << nEvents[mEntry] << " events)..." << endl;
        gRandom->SetSeed(0);

        for (unsigned i = 0; i<(unsigned)(nEvents[mEntry]); i++) {

          // Get local coordinates
          double theta = LocCoord_theta[mEntry][i];
          double phi = LocCoord_phi[mEntry][i];
          Direction dir(theta,phi);

          for (int k=0; k<opts.resample; k++) {

            // Generate new equatorial coordinates
            double rndMJD = histMJD[mEntry]->GetRandom();
            Equatorial eq = GetEquatorialFromDirection(dir, rndMJD);

            // Calculate solar dipole weighting
            double eventweight = 1.0;
            if (opts.sd)
              eventweight = solar_dipole(mjd, eq.ra, eq.dec, opts.sd2);

            // RA can change depending on time standard (anti, ext, solar, sid)
            double tmp_ra = time_standard_ra(mjd, eq, opts.method);

            // Write to map
            pointing eqDir(M_PI/2.-eq.dec, tmp_ra);
            SkyMapPtr tmp_map = BGMap[mapIdx];
            int pixelID = tmp_map->ang2pix(eqDir);
            (*tmp_map)[pixelID] += eventweight;
          }
        }
      }

      // If not on the last entry, scrambling triggered by event outside dt
      if (nevent + 1 != nEntries) {

        // Update beginning of time-scrambling window
        mjd1 += dt;
        std::cout << "new start_mjd :" << setprecision(12) << mjd1 << endl;

        // Clear storage
        for (unsigned kEntry = 0; kEntry < nMaps; kEntry++) {
          LocCoord_phi[kEntry].erase(LocCoord_phi[kEntry].begin(),
              LocCoord_phi[kEntry].end());
          LocCoord_theta[kEntry].erase(LocCoord_theta[kEntry].begin(),
              LocCoord_theta[kEntry].end());
          histMJD[kEntry]->Reset();
          nEvents[kEntry] = 0;
        }

        // Return to loop without incrementing counter to revisit event
        continue;

      }
    }

    nevent += 1;

  }


  // Finish up
  for (unsigned m=0; m<nMaps; m++) {

    SkyMap data = *(DataMap[m]);
    SkyMap bg = *(BGMap[m]);
    SkyMap local = *(LocalMap[m]);

    // Scale background map
    for (int i=0; i<bg.Npix(); i++)
      bg[i] *= alpha;

    std::cout << "Read " << nEntries << " events" << "\n"
         << "Used " << nUsedEvents[m] << " events" << endl;

    // Save BG, Data, and Local maps in one file
    arr<std::string> colname(3);
    colname[0] = "data map";
    colname[1] = "background map";
    colname[2] = "local map";

    stringstream namefits;

    namefits << opts.outdir << "/" << opts.outfile;
    if (vm.count("ebins"))
      namefits << "_" << opts.ebins[m] << "-" << opts.ebins[m+1] << "GeV";
    if (vm.count("sbins"))
      namefits << "_" << opts.sbins[m] << "to" << opts.sbins[m+1] << "s125";
    namefits << ".fits.gz";

    fitshandle fitsOut;
    fitsOut.create(namefits.str().c_str());

    fitsOut.add_comment("Maps: data, bg, local");
    //prepare_Healpix_fitsmap(fitsOut, DataMap[m],
    //    FITSUTIL<float>::DTYPE, colname);
    //    FITSUTIL<double>::DTYPE, colname);
    // Temporary workaround - Planck Data Type for double is "9"
    prepare_Healpix_fitsmap(fitsOut, data, PLANCK_FLOAT64, colname);
    fitsOut.write_column(1, data.Map());
    fitsOut.write_column(2, bg.Map());
    fitsOut.write_column(3, local.Map());
    fitsOut.close();
  }

  // Clean up
  //delete cutDST;
  for (unsigned m=0; m<nMaps; m++)
    delete histMJD[m];

  timer.Stop();
  printf("RT=%7.3f s, Cpu=%7.3f s\n",timer.RealTime(),timer.CpuTime());

}

double time_standard_ra(double mjd, Equatorial eq, std::string time_standard) {

  // Default behavior assumes sidereal time
  double ra = eq.ra;

  // Shift RA for chosen time standard
  double lst = GetGMST(mjd);
  if (time_standard == "anti") {
    double localAntiS = GetGMAST(mjd);
    ra = fmod( eq.ra - (lst + localAntiS) * M_PI/12, 2*M_PI);
  }

  else if (time_standard == "ext") {
    double localExtS = GetGMEST(mjd);
    ra = fmod( eq.ra - (lst + localExtS) * M_PI/12, 2*M_PI);
  }

  else if (time_standard == "solar") {
    double tod = ( mjd - int(mjd) ) * 24.;
    ra = fmod(eq.ra - (lst + tod) * M_PI/12., 2*M_PI);
    // Solar coordinates need an additional 180-deg flip (definition)
    ra -= M_PI;
    while (ra < 0)
      ra += 2.*M_PI;
  }

  return ra;

}



int ICenergyCut(po::variables_map vm, SimpleDST dst, photospline::splinetable<> &spline, double zenith, std::vector<float> ebins) {

  // Setup basic parameters
  double x = cos(zenith);
  double y = log10(dst.NChannels);

  // Boundary check (energy cut tables go to 0.3 in cos(zenith))
  if (x < 0.3)
    return -1;

  // Catch additional outliers
  double coords[2] = {x, y};
  int centers[spline.get_ndim()];
  if (!spline.searchcenters(coords, centers)) {
    std::cout << "Variables outside of table boundaries" << endl;
    std::cout << "x: " << x << " y: " << y << endl;
    return -1;
  }

  // Calculate reconstructed energy
  double median = spline.ndsplineeval(coords, centers, 0);
  // Make sure we're in the energy bin range
  if ((median < ebins[0]) || (median > ebins.back()))
    return -1;

  // Get energy bin
  int ebin = 0;
  while (median > ebins[ebin+1])
    ebin += 1;

  return ebin;
}

int ITs125Cut(po::variables_map vm, SimpleDST dst, std::vector<float> sbins) {

  // Get desired s125 value
  double s125 = (dst.nStations >= 5) ? dst.s125 : dst.ss125;
  double logS125 = log10(s125);

  // Make sure we're in the bin range
  if ((logS125 < sbins[0]) || (logS125 > sbins.back()))
    return -1;

  // Get s125 bin
  int sbin = 0;
  while (logS125 > sbins[sbin+1])
    sbin += 1;

  return sbin;
}


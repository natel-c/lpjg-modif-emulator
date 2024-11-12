///////////////////////////////////////////////////////////////////////////////////////
/// \file cfxinput.cpp
/// \brief extended Input module for CF conforming NetCDF files
///
/// \author Peter Anthoni
/// $Date: 2015-12-14 16:08:55 +0100 (Mon, 14 Dec 2015) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "cfxinput.h"

#ifdef HAVE_NETCDF

#include "guess.h"
#include "driver.h"
#include "weathergen.h"
#include "guessstring.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <netcdf.h>
#include <cmath>

#define MIN_GRID_DEGREES (30.0/3600.0) // 30arcsec, minimum spatial grid in deg

REGISTER_INPUT_MODULE("cfx", CFXInput)

using namespace GuessNC::CF;

namespace CFXINPUT {
XGenericSpinupData::XGenericSpinupData()
: thisyear(0) {
}

void XGenericSpinupData::get_data_from(RawData& source) {
	data = source;
	
	if (source.empty()) {
		fail("No source data given to XGenericSpinupData::get_data_from()");
	}
	
	for (size_t i = 0; i < source.size(); ++i) {
		if ((source[i].size() != 12 &&
				 source[i].size() % DAYS_PER_YEAR != 0) ||
				source[i].empty()) {
			fail("Incorrect number of timesteps per year in data given to\n"\
					 "XGenericSpinupData (got %d)", source[i].size());
		}
		
		if (i > 0 && source[i].size() != source[i-1].size()) {
			fail("Different number of timesteps in different years in\n"\
					 "data given to XGenericSpinupData::get_data_from()");
		}
	}
	
	firstyear();
}

double XGenericSpinupData::operator[](int ts) const {
	if (ts < 0 || ts >= int(data[thisyear].size())) {
		fail("Trying to access data for timestep %d during spinup\n"\
				 "(should be 0-%d)", ts, data[thisyear].size());
	}
	return data[thisyear][ts];
}

void XGenericSpinupData::nextyear() {
	thisyear = (thisyear + 1) % nbr_years();
}

void XGenericSpinupData::firstyear() {
	thisyear = 0;
}

void XGenericSpinupData::detrend_data(bool future) {
	
	std::vector<double> annual_mean(nbr_years(), 0);
	std::vector<double> year_number(nbr_years());
	
	for (size_t y = 0; y < nbr_years(); ++y) {
		for (size_t d = 0; d < data[y].size(); ++d) {
			annual_mean[y] += data[y][d];
		}
		annual_mean[y] /= data[y].size();
		year_number[y] = y;
	}
	
	double a, b;
	regress(&year_number.front(), &annual_mean.front(), nbr_years(), a, b);
	
	for (unsigned int y = 0; y < nbr_years(); ++y) {
		double anomaly;
		if(future)
			anomaly = b * (double)(y - ((int) nbr_years() - 1));
		else
			anomaly=b*(double)y;
		for (size_t d = 0; d < data[y].size(); ++d) {
			data[y][d] -= anomaly;
		}
	}
}

size_t XGenericSpinupData::nbr_years() const {
	return data.size();
}

std::string tolower(const char* str) {
	std::string result(str);
	std::transform(result.begin(), result.end(), result.begin(), ::tolower);
	return result;
}

const int SECONDS_PER_DAY = 24*60*60;

// Converts a CF standard name to one of our insolation types
// Calls fail() if the standard name is invalid
insoltype cf_standard_name_to_insoltype(const std::string& standard_name) {
	if (standard_name == "surface_downwelling_shortwave_flux_in_air" ||
	    standard_name == "surface_downwelling_shortwave_flux" ||
	    standard_name == "surface_solar_radiation") { // AB added 20191112
		return SWRAD_TS;
	}
	else if (standard_name == "surface_net_downward_shortwave_flux") {
		return NETSWRAD_TS;
	}
	else if (standard_name == "cloud_area_fraction" || standard_name == "cloud cover") {
		return SUNSHINE;
	}
	else {
		fail("Unknown insolation type: %s", standard_name.c_str());
		return SUNSHINE; // To avoid compiler warning
	}
}

// Gives the maximum allowed value for insolation, given an insolation type
// Used as an upper limit when interpolating from monthly to daily values
double max_insolation(insoltype instype) {
	if (instype == SUNSHINE) {
		return 100;
	}
	else {
		return std::numeric_limits<double>::max();
	}
}

// Checks if a DateTime is at the first day of the year
bool first_day_of_year(GuessNC::CF::DateTime dt) {
	return dt.get_month() == 1 && dt.get_day() == 1;
}

// Checks if a DateTime is in January
bool first_month_of_year(GuessNC::CF::DateTime dt) {
	return dt.get_month() == 1;
}

// Compares a Date with a GuessNC::CF::DateTime to see if the Date is on an earlier day
bool earlier_day(const Date& date, int calendar_year,
                 const GuessNC::CF::DateTime& date_time) {
	std::vector<int> d1(3),d2(3);

	d1[0] = calendar_year;
	d2[0] = date_time.get_year();

	d1[1] = date.month+1;
	d2[1] = date_time.get_month();

	d1[2] = date.dayofmonth+1;
	d2[2] = date_time.get_day();

	return d1 < d2;
}

// Compares a Date with a GuessNC::CF::DateTime to see if the Date is on a later day
// The date object must know about its calendar years (i.e. set_first_calendar_year must
// have been called)
bool later_day(const Date& date,
               const GuessNC::CF::DateTime& date_time) {
	std::vector<int> d1(3),d2(3);

	d1[0] = date.get_calendar_year();
	d2[0] = date_time.get_year();

	d1[1] = date.month+1;
	d2[1] = date_time.get_month();

	d1[2] = date.dayofmonth+1;
	d2[2] = date_time.get_day();

	return d1 > d2;
}

// Checks if the variable contains daily data
bool is_daily(const GuessNC::CF::GridcellOrderedVariable* cf_var) {

	// Check if first and second timestep is one day apart

	DateTime dt1 = cf_var->get_date_time(0);
	DateTime dt2 = cf_var->get_date_time(1);

	dt1.add_time(1, GuessNC::CF::DAYS, cf_var->get_calendar_type());

	return dt1 == dt2;
}

// Returns a DateTime in the last day for which the variable has data.
// For daily data, this is simply the day of the last timestep, for monthly data
// we need to find the last day of the last timestep's month.
GuessNC::CF::DateTime last_day_to_simulate(const GuessNC::CF::GridcellOrderedVariable* cf_var) {
	GuessNC::CF::DateTime last = cf_var->get_date_time(cf_var->get_timesteps()-1);
	if (is_daily(cf_var)) {
		return last;
	}
	else {
		// Not daily, assume monthly.
		GuessNC::CF::DateTime prev = last;
		GuessNC::CF::DateTime next = last;

		do {
			prev = next;
			next.add_time(1, GuessNC::CF::DAYS, cf_var->get_calendar_type());
		} while (next.get_month() == last.get_month());

		return prev;
	}
}

// Verifies that a CF variable with air temperature data contains what we expect
void check_temp_variable(const GuessNC::CF::GridcellOrderedVariable* cf_var) {
	if ((cf_var->get_standard_name() != "air_temperature") && (cf_var->get_long_name() != "near-surface temperature") &&
			(cf_var->get_long_name() != "2 metre temperature") && (cf_var->get_long_name() != "2m_temperature") &&
			(cf_var->get_long_name() != "temperature_at_2m")){   // added AB 20191112
		fail("Temperature variable doesn't seem to contain air CFX temperature data");
	}
	std::string units = cf_var->get_units();
	if (!((to_lower(units).compare("degrees celsius")==0) || (units.compare("K")==0) || (units.compare("C")==0))) {
		fail("Temperature variable doesn't seem to be in Kelvin, C, or degrees celsius");
	}
}

// Verifies that a CF variable with precipitation data contains what we expect
void check_prec_variable(const GuessNC::CF::GridcellOrderedVariable* cf_var) {
	if (cf_var->get_standard_name() == "precipitation_flux") {
		if (cf_var->get_units() != "kg m-2 s-1") {
			fail("Precipitation is given as flux but does not have the correct unit (kg m-2 s-1)");
		}
	}
	else if (cf_var->get_standard_name() == "precipitation_amount") {
		if (cf_var->get_units() != "kg m-2") {
			fail("Precipitation is given as amount but does not have the correct unit (kg m-2)");
		}
	}
	else if (cf_var->get_long_name() == "precipitation") {
		if (!(cf_var->get_units() == "mm/month") || (cf_var->get_units() == "mm")) { // cru ts3.24 pre has mm/month, cru_ts3.22 has mm
			fail("Precipitation does not have the correct unit (mm or mm/month)");
		}
	}
	else if (cf_var->get_long_name() == "total_precipitation") {
		if (!(cf_var->get_units() == "mm/day")) {
			fail("Precipitation does not have the correct unit (mm/day)");
		}
	}
	else {
		fail("Unrecognized precipitation type");
	}
}

// Verifies that a CF variable with insolation data contains what we expect
void check_insol_variable(const GuessNC::CF::GridcellOrderedVariable* cf_var) {
	if (cf_var->get_standard_name() != "surface_downwelling_shortwave_flux_in_air" &&
	    cf_var->get_standard_name() != "surface_downwelling_shortwave_flux" &&
	    cf_var->get_standard_name() != "surface_net_downward_shortwave_flux" &&
	    cf_var->get_standard_name() != "surface_solar_radiation" &&  // added AB 20191112
	    cf_var->get_standard_name() != "cloud_area_fraction" &&
	    cf_var->get_long_name() != "cloud cover") {
		fail("Insolation variable doesn't seem to contain insolation data");
	}

	if (cf_var->get_standard_name() == "cloud_area_fraction") {
		if (cf_var->get_units() != "1") {
			fail("Unrecognized unit for cloud cover");
		}
	} else if (cf_var->get_long_name() == "cloud cover") {
		if (cf_var->get_units() != "percentage") {
			fail("Unrecognized unit for cloud cover");
		}
	}
	else {
		std::string units = cf_var->get_units();
		if (!((units.compare("W m-2")==0) || (units.compare("W/m^2")==0))) {
			fail("Insolation variable given as radiation but unit doesn't seem to be in W m-2 or W/m^2");
		}
	}
}

// Verifies that a CF variable with wetdays data contains what we expect
void check_wetdays_variable(const GuessNC::CF::GridcellOrderedVariable* cf_var) {
	const char* wetdays_standard_name =
		"number_of_days_with_lwe_thickness_of_precipitation_amount_above_threshold";

	if (cf_var && (cf_var->get_standard_name() != wetdays_standard_name && cf_var->get_long_name() != "wet day frequency")) {
		fail("Wetdays variable should have standard name %s", wetdays_standard_name);
	}
}

// Verifies that a CF variable with pressure data contains what we expect
void check_pres_variable(const GuessNC::CF::GridcellOrderedVariable* cf_var) {
	const char* pres_standard_name = "surface_air_pressure";
	if (cf_var->get_standard_name() != pres_standard_name) {
		fail("Pressure variable should have standard name %s ",pres_standard_name);
	}
	if (cf_var->get_units() != "Pa") {
		fail("Pressure must be given in Pa!");
	}
}

// Verifies that a CF variable with specific humidity data contains what we expect
void check_specifichum_variable(const GuessNC::CF::GridcellOrderedVariable* cf_var) {
	const char* standard_name = "specific_humidity";
	if (cf_var->get_standard_name() != standard_name) {
		fail("QAir variable should have standard name %s ",standard_name);
	}
	if (cf_var->get_units() != "1") {
		fail("Specific Humidity must be dimensionless (here, '1'!");
	}
}

// Verifies that a CF variable with relative humidity data contains what we expect
void check_relhum_variable(const GuessNC::CF::GridcellOrderedVariable* cf_var) {
	const char* standard_name = "relative_humidity";
	if (cf_var->get_standard_name() != standard_name) {
		fail("Relative humidity variable should have standard name %s ",standard_name);
	}
	if (!(cf_var->get_units() == "1" || cf_var->get_units() == "" || cf_var->get_units() == "%")) {
		fail("Relative Humidity must be % or dimensionless (here, '1'!");
	}
}

// Verifies that a CF variable with wind-speed data contains what we expect
void check_wind_variable(const GuessNC::CF::GridcellOrderedVariable* cf_var) {
	const char* standard_name = "wind_speed";
	if (cf_var->get_standard_name() != standard_name) {
		fail("Wind variable should have standard name %s ",standard_name);
	}
	if (cf_var->get_units() != "m s-1") {
		fail("Wind data must be given in m s-1 !");
	}
}

// Checks if two variables contain data for the same time period
//
// Compares start and end of time series, the day numbers are only compared if
// both variables are daily.
void check_compatible_timeseries(const GuessNC::CF::GridcellOrderedVariable* var1,
                                 const GuessNC::CF::GridcellOrderedVariable* var2) {
	GuessNC::CF::DateTime start1, start2, end1, end2;

	const std::string error_message = format_string("%s and %s have incompatible timeseries",
		var1->get_variable_name().c_str(), var2->get_variable_name().c_str());

	start1 = var1->get_date_time(0);
	start2 = var2->get_date_time(0);

	end1 = var1->get_date_time(var1->get_timesteps() - 1);
	end2 = var2->get_date_time(var2->get_timesteps() - 1);

	if (start1.get_year() != start2.get_year() ||
		start1.get_month() != start2.get_month()) {
		fail(error_message.c_str());
	}

	if (end1.get_year() != end2.get_year() ||
		end1.get_month() != end2.get_month()) {
		fail(error_message.c_str());
	}

	if (is_daily(var1) && is_daily(var2)) {
		if (start1.get_day() != start2.get_day() ||
			end1.get_day() != end2.get_day()) {
			fail(error_message.c_str());
		}
	}
}

// Makes sure all variables have compatible time series
void check_compatible_timeseries(const std::vector<GuessNC::CF::GridcellOrderedVariable*> variables) {

	for (size_t i = 0; i < variables.size(); ++i) {
		for (size_t j = i + 1; j < variables.size(); ++j) {
			check_compatible_timeseries(variables[i], variables[j]);
		}
	}
}

void check_same_spatial_domains(const std::vector<GuessNC::CF::GridcellOrderedVariable*> variables) {

	for (size_t i = 1; i < variables.size(); ++i) {
		if (!variables[0]->same_spatial_domain(*variables[i])) {
			fail("%s and %s don't have the same spatial domain",
				variables[0]->get_variable_name().c_str(),
				variables[i]->get_variable_name().c_str());
		}
	}
}

bool find_closest_coord(std::vector<double> data, double pos, size_t& index, double max_dist = 360.0) {
	// max_dist setup for lon and lat to work
	bool found = false;
	//FIXME: could use sort(data) and then find lower_bound, sort requires C++17
	for (unsigned int i = 0; i < data.size(); i++) {
		double dist = fabs(data[i] - pos);
		if (dist < max_dist) {
			index = i;
			max_dist = dist;
			found = true;
		}
	}
	return found;
}

/// Retrieves dataset coordinates from (lon, lat)-coordinates
bool get_index_for_coords(double lon, double lat, size_t& x, size_t& y, std::vector<double> lons, std::vector<double> lats, double searchradius) {
	bool found = false;
	size_t lon_index = 0;
	size_t lat_index = 0;
	if (find_closest_coord(lons, lon, lon_index) && find_closest_coord(lats, lat, lat_index)) {
		x = lon_index;
		y = lat_index;
	} else {
		return false;
	}
	
	if (searchradius > 0.0) {
		if (fabs(lons[x] - lon) < searchradius && fabs(lats[y] - lat) < searchradius) {
			found = true;
		}
	} else {
		// require a direct (very very close) match
		if (fabs(lons[x]-lon)<=min(0.001,(MIN_GRID_DEGREES/2.0)) && fabs(lats[y]-lat)<=min(0.001,(MIN_GRID_DEGREES/2.0))) {
			found = true;
		}
	}

#if 0
	// previous version found a direct match or the 1st within searchradius
	// need to first look for a very close match, if not found then check within searchradius
	for (size_t d = 0; d < lons.size(); ++d) {
		if (fabs(lons[d]-lon)<=min(0.001,(MIN_GRID_DEGREES/2.0))) {
			x = d;
			found = true;
			break;
		}
	}
	if (found) {
		found = false;
		for (size_t d = 0; d < lats.size(); ++d) {
			if (fabs(lats[d]-lat)<=min(0.001,(MIN_GRID_DEGREES/2.0))) {
				y = d;
				found = true;
				break;
			}
		}
	}
	
	if (!found) {
		for (size_t d = 0; d < lons.size(); ++d) {
			if (fabs(lons[d]-lon)<=searchradius) {
				x = d;
				found = true;
				break;
			}
		}
		if (found) {
			found = false;
			for (size_t d = 0; d < lats.size(); ++d) {
				if (fabs(lats[d]-lat)<=searchradius) {
					y = d;
					found = true;
					break;
				}
			}
			
		}
	}
	
	// lets try to check if the lon+360 would work, in case the netcdf have lon 0..360
	if (!found) {
		for (size_t d = 0; d < lons.size(); ++d) {
			if (fabs(lons[d]-lon)<=searchradius || fabs(lons[d]-(lon+360.0))<=searchradius) {
				x = d;
				found = true;
				break;
			}
		}
		if (found) {
			found = false;
			for (size_t d = 0; d < lats.size(); ++d) {
				if (fabs(lats[d]-lat)<=searchradius) {
					y = d;
					found = true;
					break;
				}
			}
			
		}

	}
#endif
	return found;
}

bool cf_read_pop_dens(const GuessNC::CF::GridcellOrderedVariable* cf_popdens, std::vector<PopDensity> &pop_density) {
	int ntimes = cf_popdens->get_timesteps();
	if (ntimes == 0)
		return false;
	
	pop_density.clear();
	// lets assume we have yearly data
	PopDensity tmp_pop_density;
	for (int idx=0; idx < ntimes; idx++) {
		DateTime cf_date;
		cf_date = cf_popdens->get_date_time(idx);
		tmp_pop_density.year = cf_date.get_year();
		double density = cf_popdens->get_value(idx);
		if (std::isnan(density)||std::isinf(density)||negligible(density - (-999.0))) // GCP 2019 pop dens has -999 as missing value
			density = 0.0;
		tmp_pop_density.density = density;
		pop_density.push_back(tmp_pop_density);
	}
	
	return true;
}

// Compute relative humidity from specific humidity, temperature and pressure
double calc_relative_humidity(double temp, double specific_humidity, double pressure) {

	// qair  specific humidity, dimensionless (e.g. kg/kg) 
	// temp  temperature in degrees C
	// press pressure in Pa
	// rh    relative humidity in frac.
//	if ( pressure > 106000 || pressure < 10000 ) {
//		fail("Unit for pressure must be [Pa]: calc_relative_humidity(cfinput.cpp) year=%d Pa=%f", date.get_calendar_year(), pressure);
//	}
//	if ( temp  > 80. ) {
//		fail("Unit for temperature must be [deg C]: calc_relative_humidity(cfinput.cpp)");
//	}
//	double pres_hPa = pressure / 100.; // convert to hPa

	double pres_hPa = (pressure > 10000) ? pressure/100.0 : pressure; // might need to change from Pa to hPa

	if (temp > 80.0) // need to change from K to degC
		temp = temp - 273.15;

	// saturation water-vapour pressure following August-Roche-Magnus Formula
	double es   = 6.112 * exp(17.67 * temp/(temp + 243.5));

	// water-vapour pressure
	// derived from approximation for s = rho_w/(rho_dryAir - rho_w)    
	double e    = specific_humidity * pres_hPa / (0.378 * specific_humidity + 0.622);
	double rh   = min(max(e / es ,0.),1.) ;
	return rh;		
}

}

using namespace CFXINPUT;

CFXInput::CFXInput()
	: cf_temp(0)
	,cf_prec(0)
	,cf_insol(0)
	,cf_wetdays(0)
	,cf_min_temp(0)
	,cf_max_temp(0)
	,cf_pres(0)
	,cf_specifichum(0)
	,cf_relhum(0)
	,cf_wind(0)
	,cf_mNHxdrydep(0)
	,cf_mNOydrydep(0)
	,cf_mNHxwetdep(0)
	,cf_mNOywetdep(0)
	,cf_popdens(0)
	,ndep_timeseries("historic")
	,nodetrending(false)
	,searchradius(0)
	,fixed_spinup_climate(0)
	,fixed_ndep(0)
	,fixed_ndep_year(1850)
	,firsthistyear(-1)
	,lasthistyear(-1)
	,extend_climate(0)
	,extend_climate_nyears(NYEAR_EXTCLIM_DATA)
	,extend_climate_startyear(2099)
	,spinup_climate_nyears(NYEAR_SPINUP_DATA)
{
	// Declare instruction file parameters
	declare_parameter("ndep_timeseries", &ndep_timeseries, 10, "For Lamarque ndep: Nitrogen deposition time series to use (historic, rcp26, rcp45, rcp60 or rcp85");
	declare_parameter("nodetrending", &nodetrending, "(0) temperature data are detrended for spinup (default), (1) not detrended");
	declare_parameter("searchradius", &searchradius, 0, 100, ">0.0 use data of closest lon and lat within this search distance (degree), 0.0 (default) exact lon and lat data (within small float distance (15arcsec)) are required");
	declare_parameter("fixed_spinup_climate", &fixed_spinup_climate, 0, 1, "(0) normal run with transient climate, (1) use spinup climate for the whole run");
	declare_parameter("fixed_ndep", &fixed_ndep, 0, 3, "(0) time varying ndep (normal run, default), >0 use ndep of fixed_ndep_year (1)=for the whole run, (2)=from fixed_ndep_year onward");
	declare_parameter("fixed_ndep_year", &fixed_ndep_year, 1850, 2100, "ndep year to use for fixed_ndep>0 (def. 1850) [1850..2100]");
	
	declare_parameter("firsthistyear", &firsthistyear, -1, 10000, "First historic year after spinup");
	declare_parameter("lasthistyear", &lasthistyear, -1, 10000, "Last historic year of simulation");
	declare_parameter("extend_climate", &extend_climate, 0, 2, "Extend the climate forward (upto lasthistyear) using the last 'extend_climate_nyears' years of climate: (0) normal run w/o extending, (1)=extend run T not detrended, (2)=extend run T detrended");
	declare_parameter("extend_climate_nyears", &extend_climate_nyears, 1, 100, "Number of years to re-cycle while extending the climate (default 20) [1..100]");
	declare_parameter("extend_climate_startyear", &extend_climate_startyear, 1850, 2500, "Year from which the extend climate will start to re-cycle the climate (def. 2099, would be the last year of the repeating climate cycle) [1850..2500]");

	// NOTE: GCP only wants to use the first 20 years (e.g. 1700..1719) as spinup climate
	declare_parameter("spinup_climate_nyears", &spinup_climate_nyears, 1, 100, "Number of years to re-cycle within the LPJ-GUESS spinup (def. 30) [1..100]");
}

CFXInput::~CFXInput() {
	if (cf_temp) {
		delete cf_temp;
		cf_temp = 0;
	}
	if (cf_prec) {
		delete cf_prec;
		cf_prec = 0;
	}
	if (cf_insol) {
		delete cf_insol;
		cf_insol = 0;
	}
	if (cf_wetdays) {
		delete cf_wetdays;
		cf_wetdays = 0;
	}
	if (cf_min_temp) {
		delete cf_min_temp;
		cf_min_temp = 0;
	}
	if (cf_max_temp) {
		delete cf_max_temp;
		cf_max_temp = 0;
	}
	if (cf_pres) {
		delete cf_pres;
		cf_pres = 0;
	}
	if (cf_specifichum) {
		delete cf_specifichum;
		cf_specifichum = 0;
	}
	if (cf_relhum) {
		delete cf_relhum;
		cf_relhum = 0;
	}
	if (cf_wind) {
		delete cf_wind;
		cf_wind = 0;
	}
	if (cf_popdens) {
		delete cf_popdens;
		cf_popdens = 0;
	}
	if (cf_mNHxdrydep) {
		delete cf_mNHxdrydep;
		cf_mNHxdrydep = 0;
	}
	if (cf_mNOydrydep) {
		delete cf_mNOydrydep;
		cf_mNOydrydep = 0;
	}
	if (cf_mNHxwetdep) {
		delete cf_mNHxwetdep;
		cf_mNHxwetdep = 0;
	}
	if (cf_mNOywetdep) {
		delete cf_mNOywetdep;
		cf_mNOywetdep = 0;
	}
	
}

bool CFXInput::create_gridlist_from_cflist(xtring& file_gridlist) {

	std::ifstream ifs(file_gridlist, std::ifstream::in);

	if (!ifs.good()) {
		fail("CFXInput::init: could not open %s for input",(char*)file_gridlist);
	}

	// extract the lon lat values out of the netcdf, all netcdf need to have the same ordering!!!
	//FIXME: would be better to have a lons=cf_temp->get_lons();lats=cf_temp->get_lats(); calls, requires cfvariable.h/cpp change
	std::vector<double> lons;
	std::vector<double> lats;
	lons.empty();
	lats.empty();
	double lon,lat;
	size_t rlat=0;
	size_t max_nlons = (size_t) (360.0/MIN_GRID_DEGREES + 0.5);
	for (size_t rlon=0; rlon < max_nlons;rlon++) {
		try {
			cf_temp->get_coords_for(rlon, rlat, lon, lat);
			lons.push_back(lon);
		} catch(...) {
			// ignore em
		}
	}
	size_t rlon=0;
	for (size_t rlat=0; rlat < (size_t) (180.0/MIN_GRID_DEGREES + 0.5);rlat++) {
		try {
			cf_temp->get_coords_for(rlon, rlat, lon, lat);
			lats.push_back(lat);
		} catch(...) {
			// ignore em
		}
	}

	std::string line;
	while (getline(ifs, line)) {

		Coord c; // add new coordinate to master grid list

		// Read next record in file
		double rlat, rlon;
		int landid;
		std::string descrip;

		std::istringstream iss(line);

		bool found = false;
		if (cf_temp->is_reduced()) {
			if (iss >> landid) {
				getline(iss, descrip);

				c.landid = landid;
				// Get lon/lat for the gridcell
				cf_temp->get_coords_for(landid, c.lon, c.lat);
				c.grid_lon = c.lon;
				c.grid_lat = c.lat;
				found = true;
			}
		}
		else {
			if (iss >> rlon >> rlat) {
				// save the gridlist lon and lat, which will be used in outputs and finding other inputs (e.g. lu, soil)
				c.grid_lon = rlon;
				c.grid_lat = rlat;
				getline(iss, descrip);
				c.descrip = descrip.c_str();
				size_t ilon,ilat;
				if (get_index_for_coords(rlon, rlat, ilon, ilat, lons, lats, searchradius)) {
					dprintf("gridlist lon lat (%g,%g) index (%u,%u)\n", rlon, rlat, ilon, ilat);
					c.rlon = ilon;
					c.rlat = ilat;
					// Get netcdf lon/lat for the gridcell
					cf_temp->get_coords_for(c.rlon, c.rlat, c.lon, c.lat);
					found = true;
				} else {
					// instead of fail(), we just skip it now
					dprintf("WARNING: could not find netcdf index for (%g,%g), will be skipped\n", rlon, rlat);
					found = false;
				}
			}
		}
		if (found) {
			gridlist.push_back(c);
		}
	}
	ifs.close();

	return true;
}

void CFXInput::init() {

	// Read CO2 data from file
	co2.load_file(param["file_co2"].str);

	// Try to open the NetCDF files
	try {
		cf_temp = new GridcellOrderedVariable(param["file_temp"].str, param["variable_temp"].str);
		cf_prec = new GridcellOrderedVariable(param["file_prec"].str, param["variable_prec"].str);
		cf_insol = new GridcellOrderedVariable(param["file_insol"].str, param["variable_insol"].str);

		if (param["file_wetdays"].str != "") {
			cf_wetdays = new GridcellOrderedVariable(param["file_wetdays"].str, param["variable_wetdays"].str);
		}

		if (param["file_min_temp"].str != "") {
			cf_min_temp = new GridcellOrderedVariable(param["file_min_temp"].str, param["variable_min_temp"].str);
		}

		if (param["file_max_temp"].str != "") {
			cf_max_temp = new GridcellOrderedVariable(param["file_max_temp"].str, param["variable_max_temp"].str);
		}

		if (param["file_pres"].str != "") {
			cf_pres = new GridcellOrderedVariable(param["file_pres"].str, param["variable_pres"].str);
		}

		if (param["file_specifichum"].str != "") {
			cf_specifichum = new GridcellOrderedVariable(param["file_specifichum"].str, param["variable_specifichum"].str);
		}

		if (param["file_relhum"].str != "") {
			cf_relhum = new GridcellOrderedVariable(param["file_relhum"].str, param["variable_relhum"].str);
		}

		if (param["file_wind"].str != "") {
			cf_wind = new GridcellOrderedVariable(param["file_wind"].str, param["variable_wind"].str);
		}
		
		if (param.isparam("file_popdens")&&param["file_popdens"].str != "") {
			cf_popdens = new GridcellOrderedVariable(param["file_popdens"].str, param["variable_popdens"].str);
		}

		if (param.isparam("file_mNHxdrydep")&&param["file_mNHxdrydep"].str != "") {
			cf_mNHxdrydep = new GridcellOrderedVariable(param["file_mNHxdrydep"].str, param["variable_mNHxdrydep"].str);
		}
		
		if (param.isparam("file_mNOydrydep")&&param["file_mNOydrydep"].str != "") {
			cf_mNOydrydep = new GridcellOrderedVariable(param["file_mNOydrydep"].str, param["variable_mNOydrydep"].str);
		}
		
		if (param.isparam("file_mNHxwetdep")&&param["file_mNHxwetdep"].str != "") {
			cf_mNHxwetdep = new GridcellOrderedVariable(param["file_mNHxwetdep"].str, param["variable_mNHxwetdep"].str);
		}
		
		if (param.isparam("file_mNOywetdep")&&param["file_mNOywetdep"].str != "") {
			cf_mNOywetdep = new GridcellOrderedVariable(param["file_mNOywetdep"].str, param["variable_mNOywetdep"].str);
		}
		
	}
	catch (const std::runtime_error& e) {
		fail(e.what());
	}

	// Make sure they contain what we expect

	check_temp_variable(cf_temp);

	check_prec_variable(cf_prec);

	check_insol_variable(cf_insol);

	check_wetdays_variable(cf_wetdays);

	if (cf_min_temp) {
		check_temp_variable(cf_min_temp);
	}

	if (cf_max_temp) {
		check_temp_variable(cf_max_temp);
	}

	if (cf_pres) {
		check_pres_variable(cf_pres);
	}

	if (cf_specifichum) {
		check_specifichum_variable(cf_specifichum);
	}

	if (cf_relhum) {
		check_relhum_variable(cf_relhum);
	}

	if (cf_wind) {
		check_wind_variable(cf_wind);
	}

	check_compatible_timeseries(all_variables());

	check_same_spatial_domains(all_variables());

	if (cf_popdens) {
		if (!cf_popdens->same_spatial_domain(*cf_temp)) {
			fail("population density %s has not the same spatial domain as climate data",(char*) param["file_popdens"].str);
		};
	}
	
	if (cf_mNHxdrydep) {
		if (!cf_mNHxdrydep->same_spatial_domain(*cf_temp)) {
			fail("dry ndep %s has not the same spatial domain as climate data",(char*) param["file_mNHxdrydep"].str);
		};
	}
	
	if (cf_mNOydrydep) {
		if (!cf_mNOydrydep->same_spatial_domain(*cf_temp)) {
			fail("dry ndep %s has not the same spatial domain as climate data",(char*) param["file_mNOydrydep"].str);
		};
	}
	
	if (cf_mNHxwetdep) {
		if (!cf_mNHxwetdep->same_spatial_domain(*cf_temp)) {
			fail("wet ndep %s has not the same spatial domain as climate data",(char*) param["file_mNHxwetdep"].str);
		};
	}
	
	if (cf_mNOywetdep) {
		if (!cf_mNOywetdep->same_spatial_domain(*cf_temp)) {
			fail("wet ndep %s has not the same spatial domain as climate data",(char*) param["file_mNOywetdep"].str);
		};
	}
	

	extensive_precipitation = cf_prec->get_standard_name() == "precipitation_amount";

	// Read list of localities and store in gridlist member variable

	// Retrieve name of grid list file as read from ins file
	xtring file_gridlist=param["file_gridlist"].str;
	create_gridlist_from_cflist(file_gridlist);

	current_gridcell = gridlist.begin();

	// Open landcover files. May reduce pftlist, stlist and mtlist. Must be called before management_input->init()
	landcover_input.init();
	// Open management files
	management_input.init();
	// Open additional files
	misc_input.init();

	// state requires firsthistyear and lasthistyear to be set:
	if ((restart || save_state)) {
		if (firsthistyear<0 || lasthistyear<0) {
			fail("CFXInput: state handling requires firsthistyear and lasthistyear to be set correctly.\nEspecially if netcdf files used for the restart has a different start year than the netcdf files used for generating the state\n");
		}
	}
	
	if (firsthistyear < 0) {
		firsthistyear=cf_temp->get_date_time(0).get_year();
	}
	if (lasthistyear < 0) {
		lasthistyear = cf_temp->get_date_time(cf_temp->get_timesteps()-1).get_year();
	}
	
	//	date.set_first_calendar_year(cf_temp->get_date_time(0).get_year() - nyear_spinup);
	date.set_first_calendar_year(firsthistyear - nyear_spinup);

	soilinput.init(param["file_soildata"].str);

	if (state_year>0) {
		state_calendar_year = (state_year-1) + date.first_calendar_year;
	}

	// Set timers
	tprogress.init();
	tmute.init();

	tprogress.settimer();
	tmute.settimer(MUTESEC);
}

// helper to synchronize spinup handling
void CFXInput::spinup_data_nextyear() {
	spinup_temp.nextyear();
	spinup_prec.nextyear();
	spinup_insol.nextyear();
	if (cf_wetdays) {
		spinup_wetdays.nextyear();
	}
	if (cf_max_temp) {
		spinup_max_temp.nextyear();
	}
	if (cf_min_temp) {
		spinup_min_temp.nextyear();
	}
	if (cf_pres) {
		spinup_pres.nextyear();
	}
	if (cf_specifichum) {
		spinup_specifichum.nextyear();
	}
	if (cf_relhum) {
		spinup_relhum.nextyear();
	}
	if (cf_wind) {
		spinup_wind.nextyear();
	}
}

// helper to synchronize extended climate handling
void CFXInput::extclim_data_nextyear() {
	extclim_temp.nextyear();
	extclim_prec.nextyear();
	extclim_insol.nextyear();
	if (cf_wetdays) {
		extclim_wetdays.nextyear();
	}
	if (cf_min_temp) {
		extclim_min_temp.nextyear();
	}
	if (cf_max_temp) {
		extclim_max_temp.nextyear();
	}
	if (cf_pres) {
		extclim_pres.nextyear();
	}
	if (cf_specifichum) {
		extclim_specifichum.nextyear();
	}
	if (cf_relhum) {
		extclim_relhum.nextyear();
	}
	if (cf_wind) {
		extclim_wind.nextyear();
	}
}

bool CFXInput::getgridcell(Gridcell& gridcell) {

	double lon, lat;
	bool gridfound = false;
	// Load data for next gridcell, or if that fails, skip ahead until
	// we find one that works.
	while (current_gridcell != gridlist.end() && !gridfound) {
		
		gridfound = load_data_from_files(lon, lat);

		if (gridfound) {
			if(readdisturbance || readdisturbance_st || readelevation_st) {
				// Not all gridcells have to be included in input file
				misc_input.loaddisturbance(lon, lat);
				misc_input.loadelevation(lon, lat);
			}

		//	gridcell.climate.mean_elevation = elevation;		// Get elevation from cru_ncep
		//	if(readelevation_st)
		//		dprintf("Mean elevation = %d\n", elevation);
		}
		if (run_landcover && gridfound) {
			bool LUerror = false;
			LUerror = landcover_input.loadlandcover(current_gridcell->grid_lon, current_gridcell->grid_lat);  // we use midpoint, which is also in our gridlist_*_midpoint.txt
			if (!LUerror)
				LUerror = management_input.loadmanagement(current_gridcell->grid_lon, current_gridcell->grid_lat);
			if (LUerror) {
				//dprintf("\nError: could not find stand at (%g,%g) in landcover/management data file(s)\n", lon, lat);
				gridfound = false;
			}
		}
		if (!gridfound) {
			dprintf("\nError: could not find/read data for (%g,%g), skipping it!\n", lon, lat);
			++current_gridcell;
		}
	}
	
	if (current_gridcell == gridlist.end()) {
		// simulation finished
		return false;
	}
	
	gridcell.set_coordinates(current_gridcell->grid_lon, current_gridcell->grid_lat);
	
	// Load spinup data for all variables

	load_spinup_data(cf_temp, spinup_temp);
	load_spinup_data(cf_prec, spinup_prec);
	load_spinup_data(cf_insol, spinup_insol);

	if (cf_wetdays) {
		load_spinup_data(cf_wetdays, spinup_wetdays);
	}

	if (cf_min_temp) {
		load_spinup_data(cf_min_temp, spinup_min_temp);
	}

	if (cf_max_temp) {
		load_spinup_data(cf_max_temp, spinup_max_temp);
	}

	if (cf_pres) {
		load_spinup_data(cf_pres, spinup_pres);
	}

	if (cf_specifichum) {
		load_spinup_data(cf_specifichum, spinup_specifichum);
	}

	if (cf_relhum) {
		load_spinup_data(cf_relhum, spinup_relhum);
	}

	if (cf_wind) {
		load_spinup_data(cf_wind, spinup_wind);
	}
	
	if (!nodetrending)
		spinup_temp.detrend_data();

	// need to synchronize the spinup start climate and available climate period
	// previous handling: e.g simulation start 1700 climate start 1901, 1901 would be used in year 1200, ..., 1700, 1800, 1900
	// we need year 1920 to be used in 1200,1901 in 1201, ...!
	GuessNC::CF::DateTime first_date = cf_temp->get_date_time(0);
	int climate_first_year = first_date.get_year();
	int sync_first_year = date.first_calendar_year;
	if (restart&&((state_calendar_year+1) < climate_first_year))
		sync_first_year = state_calendar_year+1;
	int year_offset = spinup_climate_nyears-(climate_first_year - sync_first_year)%spinup_climate_nyears;

	for (int y=0;y<year_offset;y++) {
		// Move to next year in spinup dataset
		spinup_data_nextyear();
	}

	if (cf_popdens) {
		cf_read_pop_dens(cf_popdens, pop_density);
	}

	// load extend climate data for all variables
	if (extend_climate) {
		load_extclim_data(cf_temp, extclim_temp);
		load_extclim_data(cf_prec, extclim_prec);
		load_extclim_data(cf_insol, extclim_insol);

		if (extend_climate==2) {
			extclim_temp.detrend_data(true);
		}
		
		if (cf_wetdays) {
			load_extclim_data(cf_wetdays, extclim_wetdays);
		}
		
		if (cf_min_temp) {
			load_extclim_data(cf_min_temp, extclim_min_temp);
		}
		
		if (cf_max_temp) {
			load_extclim_data(cf_max_temp, extclim_max_temp);
		}

		if (cf_pres) {
			load_extclim_data(cf_pres, extclim_pres);
		}
		if (cf_specifichum) {
			load_extclim_data(cf_specifichum, extclim_specifichum);
		}
		if (cf_relhum) {
			load_extclim_data(cf_relhum, extclim_relhum);
		}
		
		if (cf_wind) {
			load_extclim_data(cf_wind, extclim_wind);
		}

		// need to synchronize the extended climate to the restart year
		// restart in history is handled in get_yearly_data
		if (restart && state_calendar_year > extend_climate_startyear) {
			// extend_climate_startyear 2012 25 yr cycle, state year 2048
			// 2013 = 1988, ..., 2049 = 1999
			// e.g. 2099 - 2099 = 0 ok, 2109 - 2099 = 10
			// e.g. 2018 - 2018 = 0 ok, 2050 - 2018 = 32 (could do a modulo)
			year_offset = state_calendar_year - extend_climate_startyear - 1; // -1 since there will be an extclim_data_nextyear later
			// we are at a restart within extclim period, need to advance to that year
			for (int y=0;y<year_offset;y++) {
				extclim_data_nextyear();
				dprintf("%d idx %d/%d cycle year %d SYNCLOOP\n",date.get_calendar_year(), extclim_temp.thisyear, extclim_temp.nbr_years(), extclim_temp.years[extclim_temp.thisyear]);
			}
		}
	}
	
	std::string insolname = cf_insol->get_standard_name();
	if (insolname == "") insolname = cf_insol->get_long_name();
	instype = cf_standard_name_to_insoltype(insolname);
	
	if (instype == NOINSOL)
		fail("Unknown insolation type: %s", insolname.c_str());
	gridcell.climate.instype = instype;
	
	dprintf("\nCommencing simulation for gridcell at (%g,%g) climate (%g,%g)\n", current_gridcell->grid_lon, current_gridcell->grid_lat, lon, lat);
	if (current_gridcell->descrip != "") {
		dprintf("Description: %s\n", current_gridcell->descrip.c_str());
	}
	
	if (param["file_ndep"].str != "") {
		double cru_lon;
		double cru_lat;
		// this tries to calculate the 0.5deg cru lon lat, might not work if that lon lat is not within the ndep data
		cru_lon = floor(current_gridcell->grid_lon * 2.0) / 2.0 + 0.25;
		cru_lat = floor(current_gridcell->grid_lat * 2.0) / 2.0 + 0.25;

		dprintf("Using Nitrogen deposition for (%3.2f,%3.2f)\n", cru_lon, cru_lat);
		// Get nitrogen deposition, using the estimated CRU coordinates
		ndep.getndep(param["file_ndep"].str, cru_lon, cru_lat,
			         Lamarque::parse_timeseries(ndep_timeseries));
	}
	
	// Setup the soil parameters
	soilinput.get_soil(current_gridcell->grid_lon, current_gridcell->grid_lat, gridcell);

	historic_timestep_temp = -1;
	historic_timestep_prec = -1;
	historic_timestep_insol = -1;
	historic_timestep_wetdays = -1;
	historic_timestep_min_temp = -1;
	historic_timestep_max_temp = -1;

	historic_timestep_pres = -1;
	historic_timestep_specifichum = -1;
	historic_timestep_relhum = -1;
	historic_timestep_wind = -1;
	

	return true;
}

bool CFXInput::load_data_from_files(double& lon, double& lat) {

	int rlon;
	int rlat;
	int landid;

	if(cf_temp->is_reduced()) {
		rlon = current_gridcell->rlon;
		rlat = current_gridcell->rlat;
		landid = current_gridcell->landid;
	}
	else {
		rlon = current_gridcell->rlon;
		rlat = current_gridcell->rlat;
	}

	// Try to load the data from the NetCDF files

	if (cf_temp->is_reduced()) {
		if (!cf_temp->load_data_for(landid) ||
		    !cf_prec->load_data_for(landid) ||
		    !cf_insol->load_data_for(landid) ||
		    (cf_wetdays && !cf_wetdays->load_data_for(landid)) ||
		    (cf_min_temp && !cf_min_temp->load_data_for(landid)) ||
		    (cf_max_temp && !cf_max_temp->load_data_for(landid)) ||
		    (cf_pres && !cf_pres->load_data_for(landid)) ||
		    (cf_specifichum && !cf_specifichum->load_data_for(landid)) ||
		    (cf_relhum && !cf_relhum->load_data_for(landid)) ||
		    (cf_wind && !cf_wind->load_data_for(landid)) ||
		    (cf_popdens && true)    || // landid not implemented for the ff. nc's, hence fail
		    (cf_mNHxdrydep && true) || 
		    (cf_mNOydrydep && true) ||
		    (cf_mNHxwetdep && true) ||
		    (cf_mNOywetdep && true)) {
			dprintf("Failed to load data for (%d) from NetCDF files, skipping.\n", landid);
			return false;
		}
	}
	else {
		if (!cf_temp->load_data_for(rlon, rlat) ||
		    !cf_prec->load_data_for(rlon, rlat) ||
		    !cf_insol->load_data_for(rlon, rlat) ||
		    (cf_wetdays && !cf_wetdays->load_data_for(rlon, rlat)) ||
		    (cf_min_temp && !cf_min_temp->load_data_for(rlon, rlat)) ||
		    (cf_max_temp && !cf_max_temp->load_data_for(rlon, rlat)) ||
		    (cf_pres && !cf_pres->load_data_for(rlon, rlat))||
		    (cf_specifichum && !cf_specifichum->load_data_for(rlon, rlat))||
		    (cf_relhum && !cf_relhum->load_data_for(rlon, rlat)) ||
		    (cf_wind && !cf_wind->load_data_for(rlon, rlat)) ||		
		    (cf_mNHxdrydep && !cf_mNHxdrydep->load_data_for(rlon, rlat)) ||
		    (cf_mNOydrydep && !cf_mNOydrydep->load_data_for(rlon, rlat)) ||
		    (cf_mNHxwetdep && !cf_mNHxwetdep->load_data_for(rlon, rlat)) ||
		    (cf_mNOywetdep && !cf_mNOywetdep->load_data_for(rlon, rlat))) {
			dprintf("Failed to load data for (%.2f,%.2f) [%d, %d] from NetCDF files, skipping.\n", current_gridcell->grid_lon, current_gridcell->grid_lat, rlon, rlat);
			return false;
		}
		
		// need to handle pop dens differently, seems GCP pop dens has some grid cells missing e.g. -105.25 73.75 pop dens is missing (-999)
		if (cf_popdens && !cf_popdens->load_data_for(rlon, rlat)) {
			// will set it to zero later
			dprintf("Failed to load pop dens for (%.2f,%.2f) [%d, %d] from NetCDF files, will use 0.0 as pop dens.\n", current_gridcell->grid_lon, current_gridcell->grid_lat, rlon, rlat);
		}

	}

	// Get lon/lat for the gridcell

	if (cf_temp->is_reduced()) {
		cf_temp->get_coords_for(landid, lon, lat);
	}
	else {
		cf_temp->get_coords_for(rlon, rlat, lon, lat);
	}

	return true;
}

void CFXInput::get_yearly_data(std::vector<double>& data,
                              const XGenericSpinupData& spinup,
                              GridcellOrderedVariable* cf_historic,
                              const XGenericSpinupData& extclim,
                              int& historic_timestep) {
	// Extract all values for this year, for one variable,
	// either from spinup dataset or historical dataset or extended climate

	int calendar_year = date.get_calendar_year();

	if (is_daily(cf_historic)) {

		data.resize(date.year_length());

		// This function is called at the first day of the year, so current_day
		// starts at Jan 1, then we step through the whole year, getting data
		// either from spinup or historical period.
		Date current_day = date;

		GuessNC::CF::DateTime first_date = cf_historic->get_date_time(0);
		GuessNC::CF::DateTime last_date = cf_historic->get_date_time(cf_historic->get_timesteps()-1);
		
		while (current_day.year == date.year) {

			// In the spinup?
			if (earlier_day(current_day, calendar_year, cf_historic->get_date_time(0)) ||
				(fixed_spinup_climate == 1)) {

				int spinup_day = current_day.day;
				// spinup object always has 365 days, deal with leap years
				if (current_day.ndaymonth[1] == 29 && current_day.month > 1) {
					--spinup_day;
				}
				data[current_day.day]  = spinup[spinup_day];
			}
			else if (calendar_year <= last_date.get_year() && !(extend_climate>0 && calendar_year > extend_climate_startyear) ) {
				// Historical period

				if (historic_timestep + 1 < cf_historic->get_timesteps()) {

					++historic_timestep;
					GuessNC::CF::DateTime dt = cf_historic->get_date_time(historic_timestep);

					if (restart) { // in case of restart from state need to advance to the restart year
						while(dt.get_year() < (state_calendar_year + 1)) {
							++historic_timestep;
							dt = cf_historic->get_date_time(historic_timestep);
						}
					}

					// Deal with calendar mismatch

					// Leap day in NetCDF variable but not in LPJ-GUESS?
					if (dt.get_month() == 2 && dt.get_day() == 29 &&
						current_day.ndaymonth[1] == 28) {
						++historic_timestep;
					}
					// Leap day in LPJ-GUESS but not in NetCDF variable?
					else if (current_day.month == 1 && current_day.dayofmonth == 28 &&
						cf_historic->get_calendar_type() == NO_LEAP) {
						--historic_timestep;
					}
				}

				if (historic_timestep < cf_historic->get_timesteps()) {
					data[current_day.day]  = cf_historic->get_value(max(0, historic_timestep));
				}
				else {
					// Past the end of the historical period, these days wont be simulated.
					data[current_day.day] = data[max(0, current_day.day-1)];
				}
			} else {
				data[current_day.day]  = extclim[current_day.day];
			}

			current_day.next();
		}
	}
	else {

		// for now, assume that data set must be monthly since it isn't daily

		data.resize(12);

		for (int m = 0; m < 12; ++m) {

			GuessNC::CF::DateTime first_date = cf_historic->get_date_time(0);
			GuessNC::CF::DateTime last_date = cf_historic->get_date_time(cf_historic->get_timesteps()-1);

			// In the spinup?
			if (calendar_year < first_date.get_year() ||
			    (calendar_year == first_date.get_year() &&
			     m+1 < first_date.get_month()) ||
			     (fixed_spinup_climate == 1)) {
				data[m] = spinup[m];
			}
			else if (calendar_year <= last_date.get_year() && !(extend_climate>0 && calendar_year > extend_climate_startyear) ) {
				// Historical period
				if (historic_timestep + 1 < cf_historic->get_timesteps()) {
					++historic_timestep;
				}
				
				GuessNC::CF::DateTime dt = cf_historic->get_date_time(historic_timestep);
				if (restart) { // in case of restart from state need to advance to the restart year
					while(dt.get_year() < (state_calendar_year +1)) {
						++historic_timestep;
						dt = cf_historic->get_date_time(historic_timestep);
					}
				}
				

				if (historic_timestep < cf_historic->get_timesteps()) {
					data[m] = cf_historic->get_value(historic_timestep);
				}
				else {
					// Past the end of the historical period, these months wont be simulated.
					data[m] = data[max(0, m-1)];
				}
			} else { // extended climate period
				data[m] = extclim[m];
			}
		}
	}
}

void CFXInput::populate_daily_array(double* daily,
                                   const XGenericSpinupData& spinup,
                                   GridcellOrderedVariable* cf_historic,
                                   const XGenericSpinupData& extclim,
                                   int& historic_timestep,
                                   double minimum,
                                   double maximum) {

	// Get the data from spinup and/or historic
	std::vector<double> data;
	get_yearly_data(data, spinup, cf_historic, extclim, historic_timestep);

	if (is_daily(cf_historic)) {
		// Simply copy from data to daily

		std::copy(data.begin(), data.end(), daily);
	}
	else {
		// for now, assume that data set must be monthly since it isn't daily

		// Interpolate from monthly to daily values

		interp_monthly_means_conserve(&data.front(), daily, minimum, maximum);
	}
}

void CFXInput::populate_daily_prec_array(long& seed) {

	// Get the data from spinup and/or historic
	std::vector<double> prec_data;
	get_yearly_data(prec_data, spinup_prec, cf_prec, extclim_prec, historic_timestep_prec);

	std::vector<double> wetdays_data;
	if (cf_wetdays) {
		get_yearly_data(wetdays_data, spinup_wetdays, cf_wetdays, extclim_wetdays, historic_timestep_wetdays);
	}

	if (is_daily(cf_prec)) {
		// Simply copy from data to daily, and if needed convert from
		// precipitation rate to precipitation amount

		for (size_t i = 0; i < prec_data.size(); ++i) {
			dprec[i] = prec_data[i];

			if (!extensive_precipitation && cf_prec->get_units() == "kg m-2 s-1") {
				dprec[i] *= SECONDS_PER_DAY;
			}
		}
	}
	else {
		// for now, assume that data set must be monthly since it isn't daily

		// TODO: change this to only run when unit is kg m-2 s-1
		// we need totals per month, so if mm/day take *date.ndaymonth[m] etc.
		// If needed convert from precipitation rate to precipitation amount
		if (!extensive_precipitation) {
			// cru ts3.x seem to have both, they are summed monthly precip
			bool mmpermonth = (cf_prec->get_units() == "mm/month") || (cf_prec->get_units() == "mm");
			for (int m = 0; m < 12; ++m) {
				if (!mmpermonth) {
  				// TODO: use the dataset's calendar type to figure out number of days in month?
				  prec_data[m] *= SECONDS_PER_DAY * date.ndaymonth[m];
				}
			}
		}

		if (cf_wetdays) {
			prdaily(&prec_data.front(), dprec, &wetdays_data.front(), seed);
		}
		else {
			interp_monthly_totals_conserve(&prec_data.front(), dprec, 0);
		}
	}
}

void CFXInput::populate_daily_arrays(Gridcell& gridcell) {
	// Extract daily values for all days in this year, either from
	// spinup dataset or historical dataset

	if ( !is_daily(cf_temp) && weathergenerator == GWGEN ) {

		int instype = cf_standard_name_to_insoltype(cf_insol->get_standard_name());
		
		// TODO IMPLEMENT cloud-frac
		if (!cf_min_temp || !cf_max_temp || !cf_wind || ( ( !cf_pres || !cf_specifichum ) && !cf_relhum ) ||
		    instype != SWRAD_TS) {
			fail("The weathergenerator GWGEN requires: \n Tmax, Tmin, (Pressure & Specific Humidity) or rel.humidity, Windspeed, and SW radiation");
		}
		
		std::vector<double> mtemp;
		get_yearly_data(mtemp, spinup_temp, cf_temp, extclim_temp, historic_timestep_temp);
		double xmtemp[12];
		double conversion = K2degC;
		if (cf_temp->get_units() != "K") conversion = 0.0; // this assumes we have already degC
		for ( int i=0; i<12; i++) {
			xmtemp[i] = mtemp[i] - conversion;
		}

		std::vector<double> mprec;
		get_yearly_data(mprec, spinup_prec, cf_prec, extclim_prec, historic_timestep_prec);
		double xmprec[12];
		double factor = 1.0;
		for ( int i=0; i<12; i++) {
			if (cf_prec->get_units() == "kg m-2 s-1")	{
				factor = SECONDS_PER_DAY * date.ndaymonth[i];
			}
			xmprec[i] = mprec[i] * factor;
		}

		std::vector<double> mwet;
		get_yearly_data(mwet, spinup_wetdays, cf_wetdays, extclim_wetdays, historic_timestep_wetdays);
		double xmwet[12];
		for ( int i=0; i<12; i++) {
			xmwet[i] = mwet[i];
		}

		std::vector<double> minsol;
		get_yearly_data(minsol, spinup_insol, cf_insol, extclim_insol, historic_timestep_insol);
		double xminsol[12];
		for ( int i=0; i<12; i++) {
			xminsol[i] = minsol[i];
		}

		std::vector<double> mtmax;
		get_yearly_data(mtmax, spinup_max_temp, cf_max_temp, extclim_max_temp, historic_timestep_max_temp);

		std::vector<double> mtmin;
		get_yearly_data(mtmin, spinup_min_temp, cf_min_temp, extclim_min_temp, historic_timestep_min_temp);

		// Record shift of mean_temperature against mean of tmin/tmax for re-adjustment
		double xmdtr[12];
		double shift[12];
		conversion = K2degC;
		if (cf_temp->get_units() != "K") conversion = 0.0; // this assumes we have already degC
		for ( int i=0; i<12; i++) {
			xmdtr[i] = 0.5 * ((mtmax[i] - conversion) - (mtmin[i] - conversion)); // Note: conversion cancels in xmdtr
			shift[i] = 0.5 * ((mtmax[i] - conversion) + (mtmin[i] - conversion)) - (mtemp[i] - conversion);
		}

		std::vector<double> mwind;
		get_yearly_data(mwind, spinup_wind, cf_wind, extclim_wind, historic_timestep_wind);
		double xmwind[12];
		for ( int i=0; i<12; i++) {
			xmwind[i] = mwind[i];
		}

		std::vector<double> mpres;
		if ( cf_pres ) {
			get_yearly_data(mpres, spinup_pres, cf_pres, extclim_pres, historic_timestep_pres);
		}

		std::vector<double> mspecifichum;
		if ( cf_specifichum ) {
			get_yearly_data(mspecifichum, spinup_specifichum, cf_specifichum, extclim_specifichum, historic_timestep_specifichum);
		}

		double xmrhum[12];
		std::vector<double> mrelhum;
		if ( cf_relhum ){ 
			get_yearly_data(mrelhum, spinup_relhum, cf_relhum, extclim_relhum, historic_timestep_relhum);
			double factor = 1.0;
			if (cf_relhum->get_units() == "%") {
				factor = 0.01;
			}
			for ( int i=0; i<12; i++) {
				xmrhum[i] = mrelhum[i] * factor;
			}
		}
		else if ( cf_pres && cf_specifichum ) {
			// compute rel. humidity if it can't be read from file
			for ( int i=0; i<12; i++) {
				xmrhum[i] = calc_relative_humidity(mtemp[i],mspecifichum[i],mpres[i]);
			}
		}

		// Use gwgen - correlated weather
		weathergen_get_met(gridcell,xmtemp,xmprec,xmwet,xminsol,xmdtr,
				   xmwind,xmrhum,dtemp,dprec,dinsol,ddtr,
				   dwind,drelhum);

		//produce tmin/tmax from daily temperature range plus shift
		int mon = 0;
		int accumday = 0;
		for (int i = 0; i < date.year_length(); ++i) {
			if ( i >= date.ndaymonth[mon]+accumday) {
				accumday += date.ndaymonth[mon];
				mon++;
			}
			dmin_temp[i] = dtemp[i] - 0.5 * ddtr[i];
			dmax_temp[i] = dtemp[i] + 0.5 * ddtr[i];
			// correct dmin and dmax against t_mean if available
			if ( cf_min_temp && cf_max_temp ) {
				dmin_temp[i] += shift[mon];
				dmax_temp[i] += shift[mon];
			}
		}
	}
	else {

		bool convT2C = (cf_temp->get_units() == "K");  // do we need to convert from K to C?
		populate_daily_array(dtemp, spinup_temp, cf_temp, extclim_temp, historic_timestep_temp, convT2C?0:-273.15);
		populate_daily_prec_array(gridcell.seed);
		populate_daily_array(dinsol, spinup_insol, cf_insol, extclim_insol, historic_timestep_insol, 0,
		                     max_insolation(instype));
		
		if (cf_min_temp) {
			populate_daily_array(dmin_temp, spinup_min_temp, cf_min_temp, extclim_min_temp, historic_timestep_min_temp, convT2C?0:-273.15);
		}
		
		if (cf_max_temp) {
			populate_daily_array(dmax_temp, spinup_max_temp, cf_max_temp, extclim_max_temp, historic_timestep_max_temp, convT2C?0:-273.15);
		}
		
		if (cf_pres) {
			populate_daily_array(dpres, spinup_pres, cf_pres, extclim_pres, historic_timestep_pres, 0);
		}
		
		if (cf_specifichum) {
			populate_daily_array(dspecifichum, spinup_specifichum, cf_specifichum, extclim_specifichum, historic_timestep_specifichum, 0);
		}
		
		if (cf_wind) {
			populate_daily_array(dwind, spinup_wind, cf_wind, extclim_wind, historic_timestep_wind, 0);
		}
		
		if (cf_relhum) {
			populate_daily_array(drelhum, spinup_relhum, cf_relhum, extclim_relhum, historic_timestep_relhum, 0);
		}
		
		// Convert to units the model expects
		for (int i = 0; i < date.year_length(); ++i) {
			if (convT2C) {
				dtemp[i] -= K2degC;
				
				if (cf_min_temp) {
					dmin_temp[i] -= K2degC;
				}
				
				if (cf_max_temp) {
					dmax_temp[i] -= K2degC;
				}
			}
			if (instype == SUNSHINE && cf_insol->get_units() == "1") {
				// Invert from cloudiness fraction (0-1) to sunshine in percent (0-100)
				dinsol[i] = (1-dinsol[i]) * 100.0;
			}
			if (instype == SUNSHINE && cf_insol->get_units() == "percent") {
				// Invert from cloudiness to sunshine in percent
				dinsol[i] = (100-dinsol[i]);
			}
						
			if ( cf_relhum && cf_relhum->get_units() == "%" ) {
				drelhum[i] = drelhum[i] * 0.01;
			}
			
			if ( (cf_pres && cf_specifichum) && !cf_relhum ) {
				// compute relative humidity for BLAZE
				drelhum[i] = calc_relative_humidity(dtemp[i], dspecifichum[i], dpres[i]);
			}
			else if ( firemodel == BLAZE && !cf_relhum ) {
				fail("BLAZE is switched on WITHOUT info on either specific humidity and pressure or relative humidity! \n" );
			}
		}
	}
	
	// Move to next year in spinup dataset
	spinup_data_nextyear();
	
	GuessNC::CF::DateTime last_date = cf_temp->get_date_time(cf_temp->get_timesteps()-1);
	if (extend_climate && date.get_calendar_year() > extend_climate_startyear) {
		if (restart||save_state) fail("sorry, cfxinput extend climate with restart from a saved state needs fixing");
		static int warnmsg=1;
		if (warnmsg) {
			dprintf("Extending the climate from %d onward by re-cycling the %d..%d years of climate.\n", extend_climate_startyear, extend_climate_startyear-extend_climate_nyears+1,extend_climate_startyear);
			warnmsg=0;
		}
		// this prints out the currently used extclim data, before the advance of the extclim for the next year (see below)
		dprintf("%d idx %d/%d cycle year %d T %7.3f P %10.6f R %10.6f\n",date.get_calendar_year(), 
			extclim_temp.thisyear, extclim_temp.nbr_years(), extclim_temp.years[extclim_temp.thisyear],
			extclim_temp[0],extclim_prec[0],extclim_insol[0]
		);
		
		if (date.get_calendar_year() > (extend_climate_startyear)) { 
			/// need to advance the extclim for next years call of populate_daily_arrays 
			extclim_data_nextyear();
		}
		// Note: extend_climate only re-cycles the climate data, not the ndep data
	}
	
	
	// Get monthly ndep values and convert to daily
	int ndep_year = date.get_calendar_year();
	if (fixed_ndep == 1) {
		ndep_year = fixed_ndep_year; // fixed ndep input year, e.g. GCP S0, S1
	}
	if (fixed_ndep == 2) { // use fixed ndep from fixed_ndep_year onward
		ndep_year = min(fixed_ndep_year,ndep_year);
	}
	if (fixed_ndep == 3) { // use fixed ndep of fixed_ndep_year for years before it
		ndep_year = max(fixed_ndep_year,ndep_year);
	}
	
	// Get monthly ndep values and convert to daily

	double mNHxdrydep[12], mNOydrydep[12];
	double mNHxwetdep[12], mNOywetdep[12];

	if (!(cf_mNHxdrydep && cf_mNOydrydep && cf_mNHxwetdep && cf_mNOywetdep)) {
		ndep.get_one_calendar_year(ndep_year,
	                               mNHxdrydep, mNOydrydep,
							       mNHxwetdep, mNOywetdep);
	} else {
		int timestep=0; // could just calculate the timestep (year-1850)*12
		GuessNC::CF::DateTime dd = cf_mNHxdrydep->get_date_time(timestep);
		GuessNC::CF::DateTime last_dd = cf_mNHxdrydep->get_date_time(cf_mNHxdrydep->get_timesteps()-1);
		if (ndep_year > last_dd.get_year()) {
			ndep_year = last_dd.get_year();
		}
		if (ndep_year > dd.get_year()) { // advance timestep to Jan in the year
			while (dd.get_year() < ndep_year) {
				timestep++;
				dd = cf_mNHxdrydep->get_date_time(timestep);
			}
		}
		
		// read 12 monthly values
		for (int m=0; m < 12; m++) {
			mNHxdrydep[m] = cf_mNHxdrydep->get_value(timestep)*86400; // convert kg m-2 s-1 -> kg m-2 d-1
			mNOydrydep[m] = cf_mNOydrydep->get_value(timestep)*86400; //
			mNHxwetdep[m] = cf_mNHxwetdep->get_value(timestep)*86400; //
			mNOywetdep[m] = cf_mNOywetdep->get_value(timestep)*86400; //
			timestep++;
		}
	}
	
	// Distribute N deposition
	distribute_ndep(mNHxdrydep, mNOydrydep,
					mNHxwetdep, mNOywetdep,
					dprec,dNH4dep,dNO3dep);
}

void CFXInput::getlandcover(Gridcell& gridcell) {

	landcover_input.getlandcover(gridcell);
	landcover_input.get_land_transitions(gridcell);
}

bool CFXInput::getclimate(Gridcell& gridcell) {

	Climate& climate = gridcell.climate;

	GuessNC::CF::DateTime last_date = last_day_to_simulate(cf_temp);

	if ((!extend_climate && (later_day(date, last_date) || date.get_calendar_year() > lasthistyear))
		|| (extend_climate>0 && date.get_calendar_year() > lasthistyear)) {
		++current_gridcell;
		return false;
	}

	climate.co2 = co2[date.get_calendar_year()];

	if (date.day == 0) {
		populate_daily_arrays(gridcell);

		// to make this work the simfire popdens needs to be disabled
		// have used some parameter.h/cpp definitions (e.g. fire_popdens_method, fixed_popdens_hist, fixed_popdens_year) for that.
		if (fire_popdens_method==2) {
			// similar to simfire code (simfire_update_pop_density)
			int cyear = date.get_calendar_year();
			if(fixed_popdens_hist>0) {
				cyear = max(fixed_popdens_year,cyear);
				if (fixed_popdens_hist==2)
					cyear = fixed_popdens_year;
			}
			
			if (pop_density.size()==0)
				fail("no population density data available, size() == 0");
			
			double popd;
			if (cyear <= pop_density[0].year) {
				popd = pop_density[0].density;
			} else {
				unsigned int idx = 0 ;
				while (pop_density[idx].year < cyear && idx < (pop_density.size()-1)) {
					idx++;
					popd=pop_density[idx].density;
				}
				if (idx >= pop_density.size()) // shouldn't really happen
					idx = pop_density.size()-1;
				if (cyear != pop_density[idx].year) { // do we need to interpolate or extrapolate?
					if (cyear > pop_density.back().year) {
						// linearly extrapolate latest growth/decline
						// FIXME: should we take a larger period for estimating the slope?
						double change_per_year = (pop_density[idx].density-pop_density[idx-1].density) / (double)((pop_density[idx].year-pop_density[idx-1].year));
						popd = pop_density[idx].density + change_per_year * (double)(cyear-pop_density[idx].year);
					} else {
						double interpf = (double)(cyear-pop_density[idx-1].year) / (double) (pop_density[idx].year - pop_density[idx-1].year);
						popd = (1. - interpf) * pop_density[idx-1].density + interpf * pop_density[idx].density;
					}
				}
			}
			
			gridcell.pop_density = max(0.,popd);
		}
	}

	climate.temp   = dtemp[date.day];
	climate.prec   = dprec[date.day];
	climate.insol  = dinsol[date.day];
	climate.relhum = drelhum[date.day];
	climate.u10    = dwind[date.day];
	climate.tmax   = dmax_temp[date.day];
	climate.tmin   = dmin_temp[date.day];
	climate.dtr    = ddtr[date.day];

	// Nitrogen deposition
	gridcell.dNH4dep = dNH4dep[date.day];
	gridcell.dNO3dep = dNO3dep[date.day];

	// bvoc
	if(ifbvoc){
		if (cf_min_temp && cf_max_temp) {
			climate.dtr = dmax_temp[date.day] - dmin_temp[date.day];
		}
		else {
			fail("When BVOC is switched on, valid paths for minimum and maximum temperature must be given.");
		}
	}

	// First day of year only ...

	if (date.day == 0) {

		// Progress report to user and update timer

		if (tmute.getprogress()>=1.0) {

			int first_historic_year = cf_temp->get_date_time(0).get_year();
			int last_historic_year = cf_temp->get_date_time(cf_temp->get_timesteps()-1).get_year();
			int historic_years = last_historic_year - first_historic_year + 1;

			int years_to_simulate = nyear_spinup + historic_years;

			int cells_done = (int)distance(gridlist.begin(), current_gridcell);

			double progress=(double)(cells_done*years_to_simulate+date.year)/
				(double)(gridlist.size()*years_to_simulate);
			tprogress.setprogress(progress);
			dprintf("%3d%% complete (%d), %s elapsed, %s remaining\n",(int)(progress*100.0),
				date.get_calendar_year(),
				tprogress.elapsed.str,tprogress.remaining.str);
			tmute.settimer(MUTESEC);
		}
	}

	return true;

}

void CFXInput::load_spinup_data(const GuessNC::CF::GridcellOrderedVariable* cf_var,
                               XGenericSpinupData& spinup_data) {

	const std::string error_message = 
		format_string("Not enough data to build spinup, at least %d years needed",
		              spinup_climate_nyears);

	XGenericSpinupData::RawData source;

	int timestep = 0;

	bool daily = is_daily(cf_var);
	// for now, assume that each data set is either daily or monthly
	bool monthly = !daily;

	// Skip the first year if data doesn't start at the beginning of the year
	while ((daily && !first_day_of_year(cf_var->get_date_time(timestep))) ||
	       (monthly && !first_month_of_year(cf_var->get_date_time(timestep)))) {
		++timestep;

		if (timestep >= cf_var->get_timesteps()) {
			fail(error_message.c_str());
		}
	}

	// Get all the values for the first NYEAR_SPINUP_DATA years,
	// and put them into source
	for (int j = 0; j < spinup_climate_nyears; ++j) {
		std::vector<double> year(daily ? GenericSpinupData::DAYS_PER_YEAR : 12);

		for (size_t i = 0; i < year.size(); ++i) {
			if (timestep < cf_var->get_timesteps()) {
				GuessNC::CF::DateTime dt = cf_var->get_date_time(timestep);

				if (daily && dt.get_month() == 2 && dt.get_day() == 29) {
					++timestep;
				}
			}

			if (timestep >= cf_var->get_timesteps()) {
				fail(error_message.c_str());
			}

			year[i] = cf_var->get_value(timestep);
			++timestep;
		}

		source.push_back(year);
	}

	spinup_data.get_data_from(source);
}

void CFXInput::load_extclim_data(const GuessNC::CF::GridcellOrderedVariable* cf_var,
															 XGenericSpinupData& extclim_data) {
	
	XGenericSpinupData::RawData source;
	
	int timestep = cf_var->get_timesteps()-1;
	
	bool daily = is_daily(cf_var);
	// for now, assume that each data set is either daily or monthly
	bool monthly = !daily;
	
	// spool backward from the end to the first year to extract
	DateTime dt = cf_var->get_date_time(timestep);
	int first_extclim_year = extend_climate_startyear - extend_climate_nyears + 1;
	while(dt.get_year() >= (first_extclim_year)) {
		dt = cf_var->get_date_time(--timestep);
	}
	dt = cf_var->get_date_time(++timestep);
	
	// Get all the values for the last extend_climate_nyears years,
	// and put them into source
	for (int i = 0; i < extend_climate_nyears; ++i) {
		std::vector<double> year(daily ? XGenericSpinupData::DAYS_PER_YEAR : 12);
		dt = cf_var->get_date_time(timestep);
		extclim_data.years.push_back(dt.get_year());
		
		for (size_t i = 0; i < year.size(); ++i) {
			year[i] = cf_var->get_value(timestep);
			++timestep;
		}
		
		source.push_back(year);
	}
	
	extclim_data.get_data_from(source);
}


namespace {

// Help function for call to remove_if below - checks if a pointer is null
bool is_null(const GuessNC::CF::GridcellOrderedVariable* ptr) {
	return ptr == 0;
}

}

std::vector<GuessNC::CF::GridcellOrderedVariable*> CFXInput::all_variables() const {
	std::vector<GuessNC::CF::GridcellOrderedVariable*> result;
	result.push_back(cf_temp);
	result.push_back(cf_prec);
	result.push_back(cf_insol);
	result.push_back(cf_wetdays);
	result.push_back(cf_min_temp);
	result.push_back(cf_max_temp);
	result.push_back(cf_pres);
	result.push_back(cf_specifichum);
	result.push_back(cf_wind);
	result.push_back(cf_relhum);

	// Get rid of null pointers
	result.erase(std::remove_if(result.begin(), result.end(), is_null),
	             result.end());

	return result;
}

#endif // HAVE_NETCDF

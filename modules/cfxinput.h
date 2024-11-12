///////////////////////////////////////////////////////////////////////////////////////
/// \file cfxiput.h
/// \brief extended Input module for CF conforming NetCDF files
///
/// \author Peter Anthoni
/// $Date: $
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_CFXInput_H
#define LPJ_GUESS_CFXInput_H

#ifdef HAVE_NETCDF

#include "soilinput.h"
#include "cruinput.h"
#include "guessnc.h"
#include <memory>
#include <limits>

namespace CFXINPUT {
	// similar to GenericSpinupData but can also detrend data used into the future
	// did not want to modify the spinup_data.h/cpp
	class XGenericSpinupData {
	public:
		static const int DAYS_PER_YEAR = 365;
		
		/// Datatype for the data, a 2D matrix of doubles
		typedef std::vector<std::vector<double> > RawData;
		
		XGenericSpinupData();
		
		/// Loads the underlying forcing data (and sets the "current" year to 0)
		void get_data_from(RawData& source);
		
		/// Gets the value for a given timestep in the "current" year
		double operator[](int ts) const;
		
		/// Goes to the next year
		void nextyear();
		
		/// Goes to the first year
		void firstyear();
		
		/// Removes trend from the original data
		void detrend_data(bool future = false);
		
		/// Returns the number of years used to construct the spinup dataset
		size_t nbr_years() const;
		
		/// for debugging: the years that are used for the spinup rawdata
		std::vector<int> years;
		
		/// The "current" year
		int thisyear;
		
	private:
		/// The forcing data which is used over and over during the spinup
		RawData data;
	};
	
	struct PopDensity {
		int year;
		double density;
	};

}

using namespace CFXINPUT;

class CFXInput : public InputModule {
public:
	CFXInput();
	
	~CFXInput();

	void init();

	/// See base class for documentation about this function's responsibilities
	bool getgridcell(Gridcell& gridcell);

	/// See base class for documentation about this function's responsibilities
	bool getclimate(Gridcell& gridcell);
	
	/// See base class for documentation about this function's responsibilities
	void getlandcover(Gridcell& gridcell);

	/// Obtains land management data for one day
	void getmanagement(Gridcell& gridcell) {
		management_input.getmanagement(gridcell, landcover_input);
	}

	/// Obtains additional environmental data that are not dynamic (e.g. elevation)
	void getmiscinput_static(Gridcell& gridcell) {
		misc_input.getmiscinput_static(gridcell);
	}

	/// Obtains additional environmental data (e.g. disturbance) for one year
	void getmiscinput_yearly(Gridcell& gridcell) {
		misc_input.getmiscinput_yearly(gridcell);
	}
	
	bool create_gridlist_from_cflist(xtring& file_gridlist);

	static const int NYEAR_SPINUP_DATA=30;
	static const int NYEAR_EXTCLIM_DATA=20;
	// to allow the spinup re-cycle period to be modified with an ins parameter, GCP wants a 20yr-repeat cycle
	int spinup_climate_nyears;

private:

	/// Land cover input module
	LandcoverInput landcover_input;

	/// Management input module
	ManagementInput management_input;

	/// Additional text data input module
	MiscInput misc_input;

	SoilInput soilinput;

	struct Coord {

		// Type for storing grid cell longitude, latitude and description text
		double grid_lon;
		double grid_lat;

		int rlon;
		int rlat;
		double lon;
		double lat;
		int landid;
		std::string descrip;
	};

	/// The grid cells to simulate
	std::vector<Coord> gridlist;

	/// The current grid cell to simulate
	std::vector<Coord>::iterator current_gridcell;

	/// search radius to use when finding netcdf data
	double searchradius;

	/// Loads data from NetCDF files for current grid cell
	/** Returns the coordinates for the current grid cell, for
	 *  the closest CRU grid cell and the soilcode for the cell.
	 *  \returns whether it was possible to load data */
	bool load_data_from_files(double& lon, double& lat);

	/// Gets the first few years of data from cf_var and puts it into spinup_data
	void load_spinup_data(const GuessNC::CF::GridcellOrderedVariable* cf_var,
	                      XGenericSpinupData& spinup_data);

	/// Gets the last few years of data from cf_var and puts it into extclim_data
	void load_extclim_data(const GuessNC::CF::GridcellOrderedVariable* cf_var,
												XGenericSpinupData& extclim_data);
	
	/// Gets data for one year, for one variable.
	/** Returns either 12 or 365/366 values (depending on LPJ-GUESS year length, not
	 *  data set year length). Gets the values from spinup and/or historic period. */
	void get_yearly_data(std::vector<double>& data,
	                     const XGenericSpinupData& spinup,
	                     GuessNC::CF::GridcellOrderedVariable* cf_historic,
	                     const XGenericSpinupData& extclim,
	                     int& historic_timestep);

	/// Fills one array of daily values with forcing data for the current year
	void populate_daily_array(double* daily,
	                          const XGenericSpinupData& spinup,
	                          GuessNC::CF::GridcellOrderedVariable* cf_historic,
	                          const XGenericSpinupData& extclim,
	                          int& historic_timestep,
	                          double minimum = -std::numeric_limits<double>::max(),
	                          double maximum = std::numeric_limits<double>::max());

	/// Same as populate_daily_array, but for precipitation which is special
	/** Uses number of wet days if available and handles extensive/intensive conversion */
	void populate_daily_prec_array(long& seed);

	/// Fills dtemp, dprec, etc. with forcing data for the current year
	void populate_daily_arrays(Gridcell& gridcell);

	/// \returns all (used) variables
	std::vector<GuessNC::CF::GridcellOrderedVariable*> all_variables() const;

//	double parse_spatial_resolution();

	/// Yearly CO2 data read from file
	/**
	 * This object is indexed with calendar years, so to get co2 value for
	 * year 1990, use co2[1990]. See documentation for GlobalCO2File for
	 * more information.
	 */
	GlobalCO2File co2;

	// The variables

	GuessNC::CF::GridcellOrderedVariable* cf_temp;

	GuessNC::CF::GridcellOrderedVariable* cf_prec;

	GuessNC::CF::GridcellOrderedVariable* cf_insol;

	GuessNC::CF::GridcellOrderedVariable* cf_wetdays;

	GuessNC::CF::GridcellOrderedVariable* cf_min_temp;

	GuessNC::CF::GridcellOrderedVariable* cf_max_temp;
	
	GuessNC::CF::GridcellOrderedVariable* cf_pres;

	GuessNC::CF::GridcellOrderedVariable* cf_specifichum;

	GuessNC::CF::GridcellOrderedVariable* cf_relhum;

	GuessNC::CF::GridcellOrderedVariable* cf_wind;

	GuessNC::CF::GridcellOrderedVariable* cf_popdens;

	// monthly mean dry nitrogen deposition (NHx kg m-2 s-1) in netCDF
	GuessNC::CF::GridcellOrderedVariable* cf_mNHxdrydep;

	// monthly mean dry nitrogen deposition (NOy, kg m-2 s-1) in netCDF
	GuessNC::CF::GridcellOrderedVariable* cf_mNOydrydep;

	// monthly mean wet nitrogen deposition (NHx, kg m-2 s-1) in netCDF
	GuessNC::CF::GridcellOrderedVariable* cf_mNHxwetdep;
	
	// monthly mean wet nitrogen deposition (NOy, kg m-2 s-1) in netCDF
	GuessNC::CF::GridcellOrderedVariable* cf_mNOywetdep;
	
	// Spinup data for each variable

	XGenericSpinupData spinup_temp;

	XGenericSpinupData spinup_prec;

	XGenericSpinupData spinup_insol;

	XGenericSpinupData spinup_wetdays;

	XGenericSpinupData spinup_min_temp;

	XGenericSpinupData spinup_max_temp;

	XGenericSpinupData spinup_pres;
	
	XGenericSpinupData spinup_specifichum;

	XGenericSpinupData spinup_relhum;
	
	XGenericSpinupData spinup_wind;
	
	// helper to synchronize spinup handling
	void spinup_data_nextyear();
	
	// extend climate bejond time in netCDF data
	/** e.g.
	extend_climate 2
	extend_climate_nyears 20
	extend_climate_startyear 2010
	lasthistyear 2050
	Extending the climate from 2010 onward by re-cycling the 1991..2010 years of climate.
	 */
	int extend_climate;
	int extend_climate_nyears;
	int extend_climate_startyear;
	
	XGenericSpinupData extclim_temp;
	XGenericSpinupData extclim_prec;
	XGenericSpinupData extclim_insol;
	XGenericSpinupData extclim_wetdays;
	XGenericSpinupData extclim_min_temp;
	XGenericSpinupData extclim_max_temp;
	XGenericSpinupData extclim_pres;
	XGenericSpinupData extclim_specifichum;
	XGenericSpinupData extclim_relhum;
	XGenericSpinupData extclim_wind;

	// helper to synchronize extended climate handling
	void extclim_data_nextyear();

	/// Temperature for current gridcell and current year (deg C)
	double dtemp[Date::MAX_YEAR_LENGTH];

	/// Precipitation for current gridcell and current year (mm/day)
	double dprec[Date::MAX_YEAR_LENGTH];

	/// Insolation for current gridcell and current year (\see instype)
	double dinsol[Date::MAX_YEAR_LENGTH];

	/// daily pressure (Pa)
	double dpres[Date::MAX_YEAR_LENGTH];

	/// daily specifichum (kg/kg)
	double dspecifichum[Date::MAX_YEAR_LENGTH];

	/// daily wind (m/s)
	double dwind[Date::MAX_YEAR_LENGTH];
	
	/// daily relative humidity (fraction)
	double drelhum[Date::MAX_YEAR_LENGTH];
	
	/// Daily N deposition for one year
	double dNH4dep[Date::MAX_YEAR_LENGTH],dNO3dep[Date::MAX_YEAR_LENGTH];

	/// Minimum temperature for current gridcell and current year (deg C)
	double dmin_temp[Date::MAX_YEAR_LENGTH];

	/// Maximum temperature for current gridcell and current year (deg C)
	double dmax_temp[Date::MAX_YEAR_LENGTH];

	/// Daily temperature range for current gridcell and current year (deg C)
	double ddtr[Date::MAX_YEAR_LENGTH];

	/// Whether the forcing data for precipitation is an extensive quantity
	/** If given as an amount (kg m-2) per timestep it is extensive, if it's
	 *  given as a mean rate (kg m-2 s-1) it is an intensive quantity */
	bool extensive_precipitation;

	// Current timestep in CF files

	int historic_timestep_temp;

	int historic_timestep_prec;

	int historic_timestep_insol;

	int historic_timestep_wetdays;

	int historic_timestep_min_temp;

	int historic_timestep_max_temp;

	int historic_timestep_pres;

	int historic_timestep_specifichum;

	int historic_timestep_relhum;

	int historic_timestep_wind;

	/// Path to CRU binary archive
//	xtring file_cru;

	/// Nitrogen deposition forcing for current gridcell
	Lamarque::NDepData ndep;

	/// Nitrogen deposition time series to use (historic,rcp26,...)
	std::string ndep_timeseries;

	// Timers for keeping track of progress through the simulation
	Timer tprogress,tmute;
	static const int MUTESEC=20; // minimum number of sec to wait between progress messages
	
	// insolation type of the radiation netcdf
	insoltype instype;
	
	// bool indicating we do not need to detrend the data
	bool nodetrending;

	// allow to keep spinup climat during the whole run
	int fixed_spinup_climate;
	// allow to keep ndep during the whole run at values of year fixed_ndep_year (def. 1860 for GCP)
	int fixed_ndep;
	int fixed_ndep_year;

	// NOTE: crop_ncep or Sam's trunk_crop_irr has a general first/lasthistyear, given as an ins parameter
	int firsthistyear;
	int lasthistyear;
	int state_calendar_year;
	

	std::vector<PopDensity> pop_density;
	
};

#endif // HAVE_NETCDF

#endif // LPJ_GUESS_CFXInput_H

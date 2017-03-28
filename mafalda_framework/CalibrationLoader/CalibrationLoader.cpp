 /*
 * 	Copyright 2011 John Idarraga, Mathieu Benoit
 *
 * 	This file is part of MAFalda.

    MAFalda is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    MAFalda is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MAFalda.  If not, see <http://www.gnu.org/licenses/>.
 */


#ifndef __CalibrationLoader_cpp
#define __CalibrationLoader_cpp

#include "CalibrationLoader.h"
#include "MAFTools.h"

using namespace MSG;

CalibrationLoader::CalibrationLoader(MediPixAlgo * alg){

	m_calibWord = 0;
	m_calibOK = false;
	m_calib = 0x0;
	m_maxWidth = 256;
	m_maxHeight = 256;

	m_callerAlgo = alg;

}

double CalibrationLoader::CalculateAndGetCalibEnergy(pair<int, int> pix, int tot) {

	// If calib is not ready.  Probably because the user didn't load the necessary files.
	// Or the CalibrationManager decided that the calibration is not good.  See output.
	if( ! CalibIsOK() ) return 0.;

	// Get the energy through calibration.
	double ec = m_calib->GetE(pix, tot);

	// As this calibrated pixel is requested fill the matrix with
	// the calibrated information in case it is needed.
	m_callerAlgo->SetCalibEnergy(pix, ec);

	return ec;
}

double CalibrationLoader::CalculateAndGetCalibEnergy(int x, int y, int tot) {
	return CalculateAndGetCalibEnergy(make_pair(x,y), tot);
}

void CalibrationLoader::SetCalibrationConfigFile_a(const char * f) {
	m_a_configFile = f; m_calibWord |= __calib_a_loaded;
	// attempt initialization of the calibration if the user finished entering the files
	if(m_calibWord == __calib_all_loaded) CalibrationInit();
}
void CalibrationLoader::SetCalibrationConfigFile_b(const char * f) {
	m_b_configFile = f; m_calibWord |= __calib_b_loaded;
	// attempt initialization of the calibration if the user finished entering the files
	if(m_calibWord == __calib_all_loaded) CalibrationInit();
}
void CalibrationLoader::SetCalibrationConfigFile_c(const char * f) {
	m_c_configFile = f; m_calibWord |= __calib_c_loaded;
	// attempt initialization of the calibration if the user finished entering the files
	if(m_calibWord == __calib_all_loaded) CalibrationInit();
}
void CalibrationLoader::SetCalibrationConfigFile_t(const char * f) {
	m_t_configFile = f; m_calibWord |= __calib_t_loaded;
	// attempt initialization of the calibration if the user finished entering the files
	if(m_calibWord == __calib_all_loaded) CalibrationInit();
}
void CalibrationLoader::SetCalibClk(double f) {
	m_calibClock = f; m_calibWord |= __calib_clock_loaded;
	if(m_calibWord == __calib_all_loaded) CalibrationInit();

}
void CalibrationLoader::CalibrationInit() {

	// If all ready to start calibration and not initilized yet
	if ( m_calibWord == __calib_all_loaded && m_calib == 0x0) {

		m_calibOK = true;

		m_calib = new MAFTools::TimePixCalibrationHandler(
				m_a_configFile.c_str(),
				m_b_configFile.c_str(),
				m_c_configFile.c_str(),
				m_t_configFile.c_str(),
				m_maxWidth, m_maxHeight);
	}

}

#endif

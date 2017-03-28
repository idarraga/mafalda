 /*
 * 	Copyright 2014 John Idarraga
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

#ifndef __DosSummaryReport_h
#define __DosSummaryReport_h

#include "TROOT.h"
#include "TString.h"

#include <vector>
#include <iostream>
#include <stdexcept>

using namespace std;

/**
 * Simple cluster container
 */
class DosCluster {

public:
	DosCluster();
	virtual ~DosCluster();

	// Sets for every cluster component
	void SetClusterPixelComponents(vector<pair<int,int > > pixels);
	void SetClusterTOTComponents(vector<int> tot);
	void SetClusterCalibEnergyComponents(vector<double> E);

	// Sets for the whole cluster
	void SetClusterDose(double E, TString name, int code = __LET_BASED_DOSE);
	void SetClusterInnerPixels(int ip){ m_clusterInnerPixels = ip; };

	enum {
		__LET_BASED_DOSE = 0
	};

private:

	// Vector of pixels coordinates composing the Cluster
	vector< pair<int, int> > m_pixels;
	// Vector of tot per pixel
	vector< int >            m_tot;
	// Vector of calib energy per pixel
	vector<double >          m_calibEnergy;

	// There is more than one way to calculate the Dose
	// In these vectors we store the dose value, the
	// name of the dose-calculation approach and a Code (see enum)
	vector< double >  m_clusterDose;
	vector< TString > m_doseTypeName;
	vector< int > m_doseTypeCode;

	int m_clusterSize;
	int m_clusterTOT;
	int m_clusterInnerPixels;

	// ROOT macro (class definition)
	// Name of the Class and Version
	ClassDef(DosCluster, 1)
};

/**
 * Serializable Dosimetry Summary Object
 */
class DosSummaryReport {

public:
	DosSummaryReport();
	// The destructor has to be virtual because of
	// the extra Methods defined by ROOT
	virtual ~DosSummaryReport();
	void PushBackDosCluster(DosCluster c){ m_clusters.push_back(c); };

	vector<DosCluster> GetDosClusterVector(){return m_clusters;};
	int GetNClusters(){return (int)m_clusters.size();};
	DosCluster GetOneCluster(int i);

	// Clear
	void Clear();

private:

	// Cluster
	vector<DosCluster> m_clusters;

	// Station Coordinates
	vector<double> m_longitude;
	vector<double> m_latitude;
	vector<double> m_altitude;

	// ROOT Macro (Class Definition)
	// Name of the Class and Version
	ClassDef(DosSummaryReport, 1)
};

#endif

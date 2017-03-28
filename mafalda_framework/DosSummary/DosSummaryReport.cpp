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

/**
 * 	These classes are a minimum dose report container
 * 	 	- DosSummaryReport: is a Dosimetry Summary Object
 *  	- DosCluster: is a cluster container
 *  Both classes are serializable into a ROOT file
 */

#include "DosSummaryReport.h"

/**
 * Serializable Dosimetry Summary Object
 */

// ROOT Macro (Class Implementation)
ClassImp(DosSummaryReport)

DosSummaryReport::DosSummaryReport(){

}

DosSummaryReport::~DosSummaryReport(){

}

DosCluster DosSummaryReport::GetOneCluster(int i) {

	try {
		DosCluster oc = m_clusters.at(i);      // vector::at throws an out-of-range
		return oc;
	} catch (out_of_range & oor) {
		cerr << "Out of Range error: " << oor.what() << ".  Return last available cluster."<< endl;
	}

	return m_clusters.at( m_clusters.size() - 1 );
}

void DosSummaryReport::Clear(){

	m_clusters.clear();

	m_latitude.clear();
	m_altitude.clear();
	m_longitude.clear();

}

/**
 * Simple cluster container
 */

// ROOT Macro (Class Implementation)
ClassImp(DosCluster)

DosCluster::DosCluster(){

}

DosCluster::~DosCluster(){

}

void DosCluster::SetClusterPixelComponents(vector<pair<int,int > > pixels){

	vector<pair<int,int > >::iterator i = pixels.begin();

	int size = 0;
	for( ; i != pixels.end() ; i++ ) {
		m_pixels.push_back( *i );
		size++;         // calculate the cluster size
	}
	m_clusterSize = size;

}

void DosCluster::SetClusterTOTComponents(vector<int> tot) {

	vector<int>::iterator i = tot.begin();
	int totalTOT = 0;
	for( ; i != tot.end() ; i++ ){
		m_tot.push_back( *i );
		totalTOT += *i; // calculate total cluster tot (cluster volume)
	}
	m_clusterTOT = totalTOT;
}

void DosCluster::SetClusterCalibEnergyComponents(vector<double> E) {

	vector<double>::iterator i = E.begin();

	for( ; i != E.end() ; i++ ){
		m_calibEnergy.push_back( *i );
	}

}

void DosCluster::SetClusterDose(double E, TString name, int code) {

	m_clusterDose.push_back( E );
	m_doseTypeName.push_back( name );
	m_doseTypeCode.push_back( code );

}


/*
 * 	Copyright 2009 John Idarraga
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

#ifndef MediPixStoreGate_cxx
#define MediPixStoreGate_cxx

#include "MPXStoreGate.h"

MediPixStoreGate::MediPixStoreGate(){

	m_nObjects = 0;

	// Set the Log service
	Log.setAlgoName("StoreGate");
	Log.OutputLevel = MSG::INFO;

}

Bool_t MediPixStoreGate::SaveObject(CandidateContainer * candidate){

	// Save pointer in the gateMap
	gateMap[candidate->GetAuthor()].push_back(candidate);

	// Keep track of the number of configuration objects
	if ( candidate->GetType() <= MPXDefs::CONF ) {
		gateMapNumberOfSpecialObjs[candidate->GetAuthor()]++;
	}

	Log << MSG::DEBUG << "Author : " << candidate->GetAuthor() << " --> " << candidate << endreq;

	return true;
}

// This function is for the user.  It delivers the numbers of objects
Int_t MediPixStoreGate::GetNObjWithAuthor(TString authorName) {

	// We deliver the number of objects available which are not of the CONFIG type
	map<TString, vector<CandidateContainer *> >::iterator i = gateMap.begin();
	vector<CandidateContainer *>::iterator icand;
	int nobj = 0;

	for( ; i != gateMap.end() ; i++ ) {
		if( (*i).first == authorName ) {
			for(icand = (*i).second.begin() ; icand != (*i).second.end() ; icand++) {
				MPXDefs::SpecialObjs type = (*icand)->GetType();
				if ( type > MPXDefs::CONF ) { // Avoid Conf objects
					nobj++;
				}
			}
			return nobj;
		}
	}

	return 0;
}

std::vector<CandidateContainer *> MediPixStoreGate::GetObjsSpecial(MPXDefs::SpecialObjs type){

	std::vector<CandidateContainer *> objs;

	std::map<TString, std::vector<CandidateContainer *> >::iterator gateItr = gateMap.begin();
	std::vector<CandidateContainer *>::iterator candItr;

	for( ; gateItr != gateMap.end() ; gateItr++) {

		candItr = (*gateItr).second.begin();

		for( ; candItr != (*gateItr).second.end() ; candItr++) {

			if ( (*candItr)->GetType() == type ) {
				objs.push_back(*candItr);
			}

		}

	}

	return objs;
}

std::vector<CandidateContainer *> MediPixStoreGate::GetObjsSpecialWithAuthor(MPXDefs::SpecialObjs type, 
		TString authorName){

	std::vector<CandidateContainer *> objs;

	std::map<TString, std::vector<CandidateContainer *> >::iterator gateItr = gateMap.begin();
	std::vector<CandidateContainer *>::iterator candItr;


	for( ; gateItr != gateMap.end() ; gateItr++)
	{
		candItr = (*gateItr).second.begin();
		for( ; candItr != (*gateItr).second.end() ; candItr++ )
			if ( (*candItr)->GetType() == type && (*candItr)->GetAuthor() == authorName)
			{
				objs.push_back(*candItr);
			}
	}

	return objs;
}

Int_t MediPixStoreGate::GetNObjWithType(MPXDefs::SpecialObjs type){

	// Between all objects in the StoreGate, regardless of their authors search for the
	//  given type.

	std::map<TString, std::vector<CandidateContainer *> >::iterator gateItr = gateMap.begin();
	std::vector<CandidateContainer *>::iterator candItr;

	Int_t total = 0;
	for( ; gateItr != gateMap.end() ; gateItr++)
	{
		candItr = (*gateItr).second.begin();
		for( ; candItr != (*gateItr).second.end() ; candItr++)
			if ( (*candItr)->GetType() == type)
				total++;
	}

	return total;
}

Int_t MediPixStoreGate::GetNObjWithTypeLessEqualThan(MPXDefs::SpecialObjs type){

	// Between all objects in the StoreGate, regardless of their authors search for the
	//  given type.

	std::map<TString, std::vector<CandidateContainer *> >::iterator gateItr = gateMap.begin();
	std::vector<CandidateContainer *>::iterator candItr;

	Int_t total = 0;
	for( ; gateItr != gateMap.end() ; gateItr++)
	{
		candItr = (*gateItr).second.begin();
		for( ; candItr != (*gateItr).second.end() ; candItr++)
			if ( (*candItr)->GetType() <= type)
				total++;
	}

	return total;
}

CandidateContainer * MediPixStoreGate::GetObjFrom(TString authorName, Int_t objIndex){

	// In case I'm asking for something that doesn't exist
	if(objIndex >= (Int_t)gateMap[authorName].size() || objIndex < 0) {
		Log << MSG::ERROR << "The object does not exist, returning pointer to NULL. The user must handle this condition." << endreq;
		return 0x0;
	}

	// The user gives the index of the work object.  I need to pass the offset of configuration
	// objects which are never erased.
	return gateMap[authorName][objIndex + gateMapNumberOfSpecialObjs[authorName] ];
}

Int_t MediPixStoreGate::CleanUpAllStoreGateExcept(MPXDefs::SpecialObjs objType){
	return CleanUpAllStoreGate(objType);
}

Int_t MediPixStoreGate::CleanUpAllStoreGate(){
	return CleanUpAllStoreGate(MPXDefs::NIL);
}

Int_t MediPixStoreGate::CleanUpAllStoreGate(MPXDefs::SpecialObjs objType){

	///////////////////////////////////////////////////////////////////
	// Some objects won't be cleared some times like --> MPXDefs::CONF
	// on an event per event basis

	///////////////////////////////
	// get number of Algos stored
	//Int_t nAlgos = (Int_t)gateMap.size();
	Int_t nObjs = 0;
	std::map<TString, std::vector<CandidateContainer *> >::iterator gateItr = gateMap.begin();
	Int_t erased = 0;

	//////////////////////////////////////
	// get number of objects for each Algo
	for(; gateItr !=  gateMap.end() ; gateItr++)
		nObjs += (Int_t)((*gateItr).second.size());

	//////////////////////////////////////
	// Delete all objects
	Int_t contItr = 0;
	for(gateItr = gateMap.begin() ; gateItr != gateMap.end() ; gateItr++)
	{
		// itr over objs on each Algo
		contItr = 0;
		//for( ; contItr < (Int_t)((*gateItr).second.size()) ; contItr++)
		while( contItr < (Int_t)((*gateItr).second.size()) )
		{

			if((*gateItr).second[contItr] != 0x0
					&& (*gateItr).second[contItr]->GetType() > objType) // > MPXDefs::CONF for example
			{
				// Erasing objects from storeGate !!!
				//  Using the member erase(..)
				//  The vector in (*gateItr).second gets shrinked by one item
				//  it means than in order to look at the next object to be erased
				//  I don't need to increment contItr.

				//(*gateItr).second.erase( (*gateItr).second.begin()+contItr );
				Log << MSG::DEBUG << "  erase ... author : " << (*gateItr).second[contItr]->GetAuthor() << endreq;
				delete (*gateItr).second[contItr];
				(*gateItr).second.erase( (*gateItr).second.begin() + contItr );


				erased++;
				//delete (*gateItr).second[contItr];
			}
			else
			{
				// if current object is not erasable, search on next item
				contItr++;
				/*
		std::cout << "[StoreGate] trying to erase a pointer to a CandidateContainer" << std::endl;
		std::cout << "            I'm not doing it ... please check you are not erasing << std::endl; 
		std::cout <<              CandidateContainer pointers yourself." << std::endl;
		std::cout << "            The StoreGate is taking care of it." << std::endl;
				 */
			}
		}

	}

	//////////////////////////////////////
	// Actually clear the store gate map.
	//  and rewind the nObjects variable.
	//gateMap.clear();

	if ( erased ) {
		Log << MSG::DEBUG << "Deleting " << erased << " objects out of " << nObjs <<endreq;
	}

	// return the actual number of Objects erased
	return erased; //nAlgos*nObjs;
}

#endif

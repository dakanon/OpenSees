/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

#ifndef PluginFramework_h
#define PluginFramework_h

// Written: Massimo Petracca 
// Created: 02/2020
// Revision: A
//
// Description: This file contains the PluginFramework manager
// class, which is a singleton that handles the loading of
// plugins

#include <PluginFrameworkAPI.h>
#include <vector>
#include <set>
#include <map>
#include <sstream>
#include <limits>
#include <iomanip>

/**
The serializer utility used to serialize heterogeneous data in a string
used to send/receive copies of plugin wrappers to/from processes 
*/
class PluginSerializer
{
private:
	std::stringstream ss;
public:
	PluginSerializer() {}
	PluginSerializer(const std::string& x)
		: ss(x) {}
	template<class T>
	PluginSerializer& write(T x) {
		ss << x << '\n';
		return *this;
	}
	PluginSerializer& write(double x) {
		ss << std::setprecision(std::numeric_limits<double>::digits10 + 1) << x << '\n';
		return *this;
	}
	PluginSerializer& write(const std::string& x) {
		ss << x << '\n';
		return *this;
	}
	template<class T>
	PluginSerializer& read(T& x) {
		std::string _x;
		std::getline(ss, _x, '\n');
		std::stringstream _ss(_x);
		_ss >> x;
		return *this;
	}
	PluginSerializer& read(std::string& x) {
		std::getline(ss, x, '\n');
		return *this;
	}
	inline std::string str()const {
		return ss.str();
	}
};

/**
Descriptor for an input argument
*/
class PluginArgumentDescriptor
{
public:
	PluginArgumentDescriptor();
	PluginArgumentDescriptor(const std::string& the_name);
	PluginArgumentDescriptor(const std::string &the_name, double def);

public:
	/// argument name
	std::string name;
	/// data for optional values
	bool is_optional;
	double default_value;
};

/**
Descriptor for a response
*/
class PluginResponseDescriptor
{
public:
	/// response name
	std::string name;
	/// all components, its length defined the size of the Vector of results
	std::vector<std::string> components;
};

class PluginLibrary; // forward declaration

/**
Descriptor for a material plugin
*/
class PluginMaterialDescriptor
{
public:
	PluginMaterialDescriptor(PluginLibrary* lib, const std::string &name, PluginMaterialProc p);
	int parseMessage(const char* m);

public:
	/// pointer to the parent library
	PluginLibrary* library;
	/// procedure name
	std::string procedure_name;
	/// pointer to the plugin material procedure
	PluginMaterialProc procedure;
	/// input arguments
	std::vector<PluginArgumentDescriptor> arguments;
	/// responses
	std::vector<PluginResponseDescriptor> responses;
private:
	bool message_parsed;
};

/**
A simple class that contains the plugin library handle, and all the functions
loaded from the library mapped to their name
*/
class PluginLibrary
{
public:
	typedef std::map<std::string, PluginMaterialDescriptor*> material_map_t;

public:
	PluginLibrary(const std::string& name, void* lh);
	~PluginLibrary();

public:
	/// the library name
	std::string library_name;
	/// the library handle
	void* library_handle;
	/// all materials loaded from this library
	material_map_t materials;
};

/**
The PluginFramework is the manager class for loading plugin libraries and procedures.
It is implemented as a singleton. Everythin is handled internally.
The only public function is getMaterialDescriptor. It may throw an exception if it cannot load
the library or the function
*/
class PluginFramework
{
private:
	typedef std::map<std::string, PluginLibrary*> plugin_map_t;

private:
	PluginFramework();
	PluginFramework(const PluginFramework&);
	PluginFramework& operator = (const PluginFramework&);

public:
	~PluginFramework();
	static PluginFramework& instance();
	///
	/// \brief getMaterialDescriptor returns a pointer to the requested PluginMaterialDescriptor
	/// \param library_name the plug-in library name
	/// \param function_name the plug-in function name name in the plug-in library
	/// \note throws a runtime_error exception in case something goes wrong while loading the library
	///
	PluginMaterialDescriptor* getMaterialDescriptor(const std::string& library_name, const std::string& function_name);
	///
	/// \brief makeMaterialData creates a new material data
	/// and initializes its fields to default values.
	/// \param p the plug-in function pointer
	/// \param tag the material tag
	///
	PluginMaterialData* makeMaterialData(PluginMaterialProc p, int tag);
	///
	/// \brief allocateData allocates the fields in d based on the material type specified in d
	/// during the first call to the material procedure with job = PF_MAT_INITIALIZE
	/// \param d the data
	///
	int allocateData(PluginMaterialData* d);
	///
	/// \brief releaseData destroys the fields in d previously allocated by a call to
	/// \ref allocateData. This method should be called after the call to the material
	/// procedure with job = PF_MAT_FINALIZE
	/// \param d the data
	///
	int releaseData(PluginMaterialData* d);
	///
	/// \brief parseTclCommand parses the TCL command according to the arguments in
	/// the descriptor, and places the obtained values in the data param vector
	///
	int parseTclCommand(const PluginMaterialDescriptor* descriptor, PluginMaterialData* d);

private:
	/// all the loaded plugins
	plugin_map_t plugins;
};

#endif // PluginFramework_h


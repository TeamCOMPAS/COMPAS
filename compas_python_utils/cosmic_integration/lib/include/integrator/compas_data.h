#ifndef _INTEGRATOR_COMPAS_DATA_READER_H_
#define _INTEGRATOR_COMPAS_DATA_READER_H_
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <H5Cpp.h>

namespace integrator{ namespace compas_data_reader{
    // This class loads data from a COMPAS HDF5 file
    class CompasData {

        private:
            double* load_data(H5::H5File file, std::string dataset_name);
            void SetMasks();
            void FindStarFormingMassPerBinary();

        public:
            CompasDataReader(std::string filename);
            ~CompasDataReader();

            double* m1_zams;
            double* m2_zams;
            double* m1;
            double* m2;
            
            int n_systems;
            int n_binaries;

    };

    void CompasData::CompasDataReader(std::string filename)
    {
        // Read the COMPAS HDF5 file and create a CompasData object
        // This function is called by the constructor

        // Open the file
        H5::H5File file(filename, H5F_ACC_RDONLY);

        // Read in data
        m1_zams = load_data(file, "BSE_System_Parameters/Mass@ZAMS(1)");
        m2_zams = load_data(file, "BSE_System_Parameters/Mass@ZAMS(2)");
        m1 = load_data(file, "BSE_Double_Compact_Objects/Mass(1)");
        m2 = load_data(file, "BSE_Double_Compact_Objects/Mass(2)");
        formation_time = load_data(file, "BSE_Double_Compact_Objects/Time");
        coalescence_time = load_data(file, "BSE_Double_Compact_Objects/Coalescence_Time");
        dco_seeds = load_data(file, "BSE_Double_Compact_Objects/SEED");
        // ...
        // Close the file
    }






}}
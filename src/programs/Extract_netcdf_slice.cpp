#include "Mesh.h"
#include <cstdlib>
int main(int argc, char *argv[])
{

    if (argc < 5)
    {
        printf("Usage: %s input_model ouput_slice slice_direction position\n", argv[0]);
        return 1;
    }
    char *filename = argv[1];
    char *out_filename = argv[2];
    cout << "filename: " << filename << endl;
    cout << "out_filename: " << out_filename << endl;
    string axis = argv[3];
    double r0, lat0, lon0;
    unsigned int ind_r, ind_lat, ind_lon;
    if (axis == "lat")
    {
        // ind_lat = atoi(argv[4]);
        string temp(argv[4]);
        istringstream iss(temp);
        iss >> lat0;
        ind_lon = 0;
        ind_r = 0;
        cout<<"latitude: "<<lat0<<endl;
        // cout << "ind_lat " << ind_lat << endl;
    }
    else if (axis == "lon")
    {
        string temp(argv[4]);
        istringstream iss(temp);
        iss >> lon0;

        ind_lat = 0;
        ind_r = 0;
        cout<<"longitude: "<<lon0<<endl;
        // cout << "ind_lon " << ind_lon << endl;
    }
    else if (axis == "r")
    {
        string temp(argv[4]);
        istringstream iss(temp);
        iss >> r0;
        ind_lat = 0;
        ind_lon = 0;
        cout<<"r: "<<r0<<endl;
        // cout << "ind_r " << ind_r << endl;
    }
    else
    {
        cout << "The third parameter must be: \"lat\" or \"lon\" or \"r\"" << endl;
    }

    vector<vector<double>> slice;

    const char *LAT_NAME = "latitude";
    const char *LON_NAME = "longitude";
    const char *DEN_NAME = (argc==6)?(argv[5]):("density");
    cout<<"variable name: "<<DEN_NAME<<endl;
    const char *R_NAME = "radius";

    string UNITS = "units";
    string POSITIVE = "positive";
    string DEGREES_EAST = "degrees_east";
    string DEGREES_NORTH = "degrees_north";
    string KILOMETER = "kilometer";
    string UP = "up";
    // For the units attributes.
    string DEN_UNITS = "kg_per_cubic_meter";
    string LAT_UNITS = "degrees_north";
    string LON_UNITS = "degrees_east";
    string BOUNDS = "bounds";
    string LATBND_NAME = "latbnd";
    string LONBND_NAME = "lonbnd";
    string RBND_NAME = "rbnd";
    double *lats;
    double *lons;
    double *rs;
    unsigned int NLAT;
    unsigned int NLON;
    unsigned int NR;
    double min_lat, max_lat, min_lon, max_lon, min_r, max_r;
    double *range_r, *range_lat, *range_lon;
    try
    {
        NcFile dataFile(filename, NcFile::read);
        NcVar data = dataFile.getVar(DEN_NAME);
        if (data.isNull())
        {
            cout << "Data is null." << endl;
            return NC_ERR;
        }
        // cout << "A0" << endl;
        vector<NcDim> dimVector_in = data.getDims();
        NcDim latDim_in = dimVector_in[0];
        NcDim lonDim_in = dimVector_in[1];
        NcDim rDim_in = dimVector_in[2];
        NLAT = latDim_in.getSize();
        NLON = lonDim_in.getSize();
        NR = rDim_in.getSize();

        cout << "NLAT=" << NLAT << endl;
        cout << "NLON=" << NLON << endl;
        cout << "NR=" << NR << endl;
        cout << "ind_lat=" << ind_lat << endl;
        cout << "ind_lon=" << ind_lon << endl;
        cout << "ind_r=" << ind_r << endl;

        lats = new double[NLAT];
        lons = new double[NLON];
        rs = new double[NR];

        NcVar latbnd = dataFile.getVar("latbnd");
        NcVar lonbnd = dataFile.getVar("lonbnd");
        NcVar rbnd = dataFile.getVar("rbnd");

        vector<size_t> startp_bnd, countp_bnd;
        startp_bnd.push_back(0);
        startp_bnd.push_back(0);
        countp_bnd.push_back(1);
        countp_bnd.push_back(1);

        double a, b;
        for (int i = 0; i < NLAT; i++)
        {
            startp_bnd[0] = i;
            startp_bnd[1] = 0;
            latbnd.getVar(startp_bnd, countp_bnd, &a);

            startp_bnd[0] = i;
            startp_bnd[1] = 1;
            latbnd.getVar(startp_bnd, countp_bnd, &b);

            lats[NLAT - 1 - i] = 0.5 * (a + b);
            if (i == 0)
            {
                max_lat = a;
            }
            if (i == NLAT - 1)
            {
                min_lat = b;
                // lats[0] = b;
            }
        }
        range_lat = new double[2];
        range_lat[0] = min_lat;
        range_lat[1] = max_lat;
        // for (int i = 0; i < NLAT; i++)
        // {
        //     cout << "lats" << i << " " << lats[i] << endl;
        // }
        // cout << "max_lat" << max_lat << ", min_lat" << min_lat << endl;
        // startp_bnd[1] = 1;
        // latbnd.getVar(startp_bnd, countp_bnd, &a);
        // lats[NLAT] = a;

        for (int i = 0; i < NLON; i++)
        {
            startp_bnd[0] = i;
            startp_bnd[1] = 0;
            lonbnd.getVar(startp_bnd, countp_bnd, &a);

            startp_bnd[0] = i;
            startp_bnd[1] = 1;
            lonbnd.getVar(startp_bnd, countp_bnd, &b);
            lons[i] = 0.5 * (a + b);
            if (i == 0)
            {
                min_lon = a;
            }
            if (i == NLON - 1)
            {
                max_lon = b;
            }
        }
        range_lon = new double[2];
        range_lon[0] = min_lon;
        range_lon[1] = max_lon;

        // for (int i = 0; i < NLON; i++)
        // {
        //     cout << "lons" << i << " " << lons[i] << endl;
        // }
        // startp_bnd[1] = 1;
        // lonbnd.getVar(startp_bnd, countp_bnd, &a);
        // lons[NLON] = a;
        // cout << "A0" << endl;

        for (int i = 0; i < NR; i++)
        {
            startp_bnd[0] = i;
            startp_bnd[1] = 0;
            rbnd.getVar(startp_bnd, countp_bnd, &a);
            a=a/1000.0;

            startp_bnd[0] = i;
            startp_bnd[1] = 1;
            rbnd.getVar(startp_bnd, countp_bnd, &b);
            b=b/1000.0;
            rs[i] = 0.5 * (a + b);
            if (i == 0)
            {
                min_r = a;
            }
            if (i == NR - 1)
            {
                max_r = b;
            }
        }
        range_r = new double[2];
        range_r[0] = min_r;
        range_r[1] = max_r;
        // for (int i = 0; i < NR; i++)
        // {
        //     cout << "rs" << i << " " << rs[i] << endl;
        // }
        // cout << "min_r " << min_r << ", max_r" << max_r << endl;
        // cout << endl;
        // startp_bnd[1] = 1;
        // rbnd.getVar(startp_bnd, countp_bnd, &a);
        // rs[NR] = a;
        // NcVar latVar = dataFile.getVar(LAT_NAME);
        // NcVar lonVar = dataFile.getVar(LON_NAME);
        // NcVar rVar = dataFile.getVar(R_NAME);
        // cout << "A" << endl;
        if (axis == "lat")
        {
            cout<<"slice along latitude"<<endl;
            double a, b;
            for (int i_lat = 0; i_lat < NLAT; i_lat++)
            {
                startp_bnd[0] = i_lat;
                startp_bnd[1] = 0;
                latbnd.getVar(startp_bnd, countp_bnd, &a);

                startp_bnd[0] = i_lat;
                startp_bnd[1] = 1;
                latbnd.getVar(startp_bnd, countp_bnd, &b);
                // cout<<a<<","<<b<<endl;
                if(a>b){
                    double temp=a;
                    a=b;
                    b=temp;
                }
                if ((lat0 > a || std::abs(lat0 - a) < 1e-10) && (lat0 < b || std::abs(lat0 - b) < 1e-10))
                {
                    ind_lat = i_lat;
                    cout<<"ind_lat: "<<ind_lat<<endl;
                    break;
                }
            }
            assert(ind_lat < NLAT);
            slice.resize(NR);
            for (int i = 0; i < NR; i++)
            {
                slice[i].resize(NLON);
            }
            vector<size_t> startp, countp;
            startp.push_back(ind_lat);
            startp.push_back(0);
            startp.push_back(0);

            countp.push_back(1);
            countp.push_back(1);
            countp.push_back(1);
            for (size_t i = 0; i < NLON; i++)
            {
                for (size_t j = 0; j < NR; j++)
                {
                    startp[1] = i;
                    startp[2] = j;
                    double a;
                    data.getVar(startp, countp, &a);
                    slice[j][i] = a;
                }
            }
            try
            {
                NcFile test(out_filename, NcFile::replace);
                test.putAtt("Conventions", "CF-1.5");
                // test.putAtt("node_offset", "1"); //0 for gridline node registration (default), 1 for pixel registration
                test.putAtt("node_offset", ncInt, 1);
                // NcDim latDim2 = test.addDim(LAT_NAME, NLAT);
                // NcVar latVar = test.addVar(LAT_NAME, ncDouble, latDim2);
                // latVar.putAtt("long_name", "Latitude");
                // latVar.putAtt("units", DEGREES_NORTH);
                // latVar.putVar(lats);

                NcDim rDim2 = test.addDim(R_NAME, NR);
                NcVar rVar = test.addVar(R_NAME, ncDouble, rDim2);
                rVar.putAtt("long_name", "Radius");
                rVar.putAtt("units", KILOMETER);
                rVar.putAtt("actual_range", ncDouble, 2, range_r);
                rVar.putVar(rs);

                NcDim lonDim2 = test.addDim(LON_NAME, NLON);
                NcVar lonVar = test.addVar(LON_NAME, ncDouble, lonDim2);
                lonVar.putAtt("long_name", "Longitude");
                lonVar.putAtt("units", DEGREES_EAST);
                lonVar.putAtt("actual_range", ncDouble, 2, range_lon);
                lonVar.putVar(lons);

                vector<NcDim> dimVector2;
                dimVector2.push_back(rDim2);
                dimVector2.push_back(lonDim2);
                NcVar denVar = test.addVar(DEN_NAME, ncDouble, dimVector2);

                vector<size_t>
                    startp2,
                    countp2;
                startp2.push_back(0);
                startp2.push_back(0);

                countp2.push_back(1);
                countp2.push_back(1);

                for (size_t i = 0; i < NR; i++)
                {
                    for (size_t j = 0; j < NLON; j++)
                    {

                        startp2[0] = i;
                        startp2[1] = j;
                        double a = slice[i][j];
                        denVar.putVar(startp2, countp2, &a);
                    }
                }
                delete[] lats;
                lats = NULL;

                delete[] lons;
                lons = NULL;

                delete[] rs;
                rs = NULL;
                              
                delete[] range_r;
                range_r=NULL;
                delete[] range_lat;
                range_lat=NULL;
                delete[] range_lon;
                range_lon=NULL;                
                cout << "The slice has been written to NetCDF file: " << out_filename << endl;
            }
            catch (NcException &e)
            {
                e.what();
                cout << "here3" << endl;
                return NC_ERR;
            }
        }
        else if (axis == "lon")
        {
            double a, b;
            for (int i_lon = 0; i_lon < NLON; i_lon++)
            {
                startp_bnd[0] = i_lon;
                startp_bnd[1] = 0;
                lonbnd.getVar(startp_bnd, countp_bnd, &a);

                startp_bnd[0] = i_lon;
                startp_bnd[1] = 1;
                lonbnd.getVar(startp_bnd, countp_bnd, &b);
                if ((lon0 > a || std::abs(lon0 - a) < 1e-10) && (lon0 < b || std::abs(lon0 - b) < 1e-10))
                {
                    ind_lon = i_lon;
                    cout<<"ind_lon: "<<ind_lon<<endl;
                    break;
                }
            }
            assert(ind_lon < NLON);
            slice.resize(NR);
            for (int i = 0; i < NR; i++)
            {
                slice[i].resize(NLAT);
            }
            vector<size_t> startp, countp;
            startp.push_back(0);
            startp.push_back(ind_lon);
            startp.push_back(0);

            countp.push_back(1);
            countp.push_back(1);
            countp.push_back(1);
            for (size_t i = 0; i < NLAT; i++)
            {
                for (size_t j = 0; j < NR; j++)
                {
                    startp[0] = i;
                    startp[2] = j;
                    double a;
                    data.getVar(startp, countp, &a);
                    slice[j][NLAT - 1 - i] = a;
                }
            }
            try
            {
                NcFile test(out_filename, NcFile::replace);
                test.putAtt("Conventions", "CF-1.5");
                test.putAtt("node_offset", ncInt, 1); //0 for gridline node registration (default), 1 for pixel registration
                NcDim latDim2 = test.addDim(LAT_NAME, NLAT);
                NcVar latVar = test.addVar(LAT_NAME, ncDouble, latDim2);
                latVar.putAtt("long_name", "Latitude");
                latVar.putAtt("units", DEGREES_NORTH);
                latVar.putAtt("actual_range", ncDouble, 2, range_lat);
                latVar.putVar(lats);

                NcDim rDim2 = test.addDim(R_NAME, NR);
                NcVar rVar = test.addVar(R_NAME, ncDouble, rDim2);
                rVar.putAtt("long_name", "Radius");
                rVar.putAtt("units", KILOMETER);
                rVar.putAtt("actual_range", ncDouble, 2, range_r);
                rVar.putVar(rs);
                // NcDim lonDim2 = test.addDim(LON_NAME, NLON);
                // NcVar lonVar = test.addVar(LON_NAME, ncDouble, lonDim2);
                // lonVar.putAtt("long_name", "Longitude");
                // lonVar.putAtt("units", DEGREES_EAST);
                // lonVar.putVar(lons);

                vector<NcDim> dimVector2;
                dimVector2.push_back(rDim2);
                dimVector2.push_back(latDim2);
                NcVar denVar = test.addVar(DEN_NAME, ncDouble, dimVector2);

                vector<size_t>
                    startp2,
                    countp2;
                startp2.push_back(0);
                startp2.push_back(0);

                countp2.push_back(1);
                countp2.push_back(1);

                for (size_t i = 0; i < NR; i++)
                {
                    for (size_t j = 0; j < NLAT; j++)
                    {

                        startp2[0] = i;
                        startp2[1] = j;
                        double a = slice[i][j];
                        denVar.putVar(startp2, countp2, &a);
                    }
                }
                delete[] lats;
                lats = NULL;

                delete[] lons;
                lons = NULL;

                delete[] rs;
                rs = NULL;
                
                delete[] range_r;
                range_r=NULL;
                delete[] range_lat;
                range_lat=NULL;
                delete[] range_lon;
                range_lon=NULL;
                cout << "The slice has been written to NetCDF file: " << out_filename << endl;
            }
            catch (NcException &e)
            {
                e.what();
                cout << "here3" << endl;
                return NC_ERR;
            }
        }
        else if (axis == "r")
        {
            double a, b;
            for (int i_r = 0; i_r < NR; i_r++)
            {
                startp_bnd[0] = i_r;
                startp_bnd[1] = 0;
                rbnd.getVar(startp_bnd, countp_bnd, &a);

                startp_bnd[0] = i_r;
                startp_bnd[1] = 1;
                rbnd.getVar(startp_bnd, countp_bnd, &b);
                if ((r0 > a || std::abs(r0 - a) < 1e-10) && (r0 < b || std::abs(r0 - b) < 1e-10))
                {
                    ind_r=i_r;
                    cout<<"ind_r: "<<ind_r<<endl;
                }
            }

            assert(ind_r < NR);
            slice.resize(NLAT);
            for (int i = 0; i < NLAT; i++)
            {
                slice[i].resize(NLON);
            }
            vector<size_t> startp, countp;
            startp.push_back(0);
            startp.push_back(0);
            startp.push_back(ind_r);

            countp.push_back(1);
            countp.push_back(1);
            countp.push_back(1);
            for (size_t i = 0; i < NLAT; i++)
            {
                for (size_t j = 0; j < NLON; j++)
                {
                    startp[0] = i;
                    startp[1] = j;
                    double a;
                    data.getVar(startp, countp, &a);
                    slice[NLAT - 1 - i][j] = a;
                }
            }
            try
            {
                NcFile test(out_filename, NcFile::replace);
                test.putAtt("Conventions", "CF-1.5");
                test.putAtt("node_offset", ncInt, 1); //0 for gridline node registration (default), 1 for pixel registration
                NcDim latDim2 = test.addDim(LAT_NAME, NLAT);
                NcVar latVar = test.addVar(LAT_NAME, ncDouble, latDim2);
                latVar.putAtt("long_name", "Latitude");
                latVar.putAtt("units", DEGREES_NORTH);

                latVar.putAtt("actual_range", ncDouble, 2, range_lat);
                // latVar.putAtt("minimum", ncDouble, min_lat);
                // latVar.putAtt("maximum", ncDouble, max_lat);
                latVar.putVar(lats);

                NcDim lonDim2 = test.addDim(LON_NAME, NLON);
                NcVar lonVar = test.addVar(LON_NAME, ncDouble, lonDim2);
                lonVar.putAtt("long_name", "Longitude");
                lonVar.putAtt("units", DEGREES_EAST);
                lonVar.putAtt("actual_range", ncDouble, 2, range_lon);

                // lonVar.putAtt("minimum", ncDouble, min_lon);
                // lonVar.putAtt("maximum", ncDouble, max_lon);
                // lonVar.putAtt("range", ncDouble, 2, range_lon);
                lonVar.putVar(lons);

                vector<NcDim> dimVector2;
                dimVector2.push_back(latDim2);
                dimVector2.push_back(lonDim2);
                NcVar denVar = test.addVar(DEN_NAME, ncDouble, dimVector2);
                vector<size_t>
                    startp2,
                    countp2;
                startp2.push_back(0);
                startp2.push_back(0);

                countp2.push_back(1);
                countp2.push_back(1);

                for (size_t i = 0; i < NLAT; i++)
                {
                    for (size_t j = 0; j < NLON; j++)
                    {

                        startp2[0] = i;
                        startp2[1] = j;
                        double a = slice[i][j];
                        denVar.putVar(startp2, countp2, &a);
                    }
                }
                delete[] lats;
                lats = NULL;

                delete[] lons;
                lons = NULL;

                delete[] rs;
                rs = NULL;
                
                delete[] range_r;
                range_r=NULL;
                delete[] range_lat;
                range_lat=NULL;
                delete[] range_lon;
                range_lon=NULL;
                cout << "The slice has been written to NetCDF file: " << out_filename << endl;
            }
            catch (NcException &e)
            {
                e.what();
                cout << "here3" << endl;
                return NC_ERR;
            }
        }
    }
    catch (NcException &e)
    {
        e.what();
        cout << "here2" << endl;
        return NC_ERR;
    }

    return 0;
}

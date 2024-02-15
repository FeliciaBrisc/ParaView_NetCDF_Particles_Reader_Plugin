""" 
If you didn't already set the Python/System path for the netcdf4-python package
then set it up here 
#import sys  
#sys.path.append('/...') 
"""


import netCDF4
import vtk  
from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
 


# The paraview.util.vtkAlgorithm module provides VTKPythonAlgorithmBase, the base class
# for all python-based vtkAlgorithm subclasses in VTK and decorators used to
# 'register' the algorithm with ParaView along with information about the GUI.
from paraview.util.vtkAlgorithm import *


"""
Code by Felicia Brisc (CEN University of Hamburg), distributed under a MIT License
A Python filter for ParaView (www.paraview.org). Reads particles trajectories files in NetCDF format
Version 1.0 
The reader requires the external module netcdf4-python
The examples made available by Kitware at the link below have provided the starting point for this reader
https://gitlab.kitware.com/paraview/paraview/blob/master/Examples/Plugins/PythonAlgorithm/PythonAlgorithmExamples.py
"""


def createModifiedCallback(anobject):
    import weakref
    weakref_obj = weakref.ref(anobject)
    anobject = None
    def _markmodified(*args, **kwars):
        o = weakref_obj()
        if o is not None:
            o.Modified()
    return _markmodified



# Decorators used to add the reader 
# @smproxy.source(name="PythonNetCDFProfilesReader", label="Python-based NetCDF Profiles Reader")
# @smhint.xml("""<ReaderFactory extensions="nc" file_description="Numpy NetCDF Profiles files" />""")
@smproxy.reader(name="PythonNetCDFParticlesReader", label="Python-based NetCDF Particles Reader",
                extensions="nc", file_description="NetCDF Particles files")



class PythonNetCDFProfilesReader(VTKPythonAlgorithmBase):
    """A reader that reads a NetCDF vertical profile file. The NetCDF file needs to  a "time" dimension, so 
    the data is treated as a temporal dataset"""

    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=0, nOutputPorts=1, outputType='vtkPolyData')
        self.filename = None
        self.nc_dataset = None

        self.newPoints = vtk.vtkPoints()

        #this will form the base for the ParaView TemporalParticlesToPathlines filter
        self.pdo = vtk.vtkPolyData()
        self.polyIndex = 0
 
        #the mandatory dimensions - will be put separately from the other (generic) variables  
        self.lat = None
        self.lon = None        
        self.time = None
        
        self.numPoints = 0

        self.sphericalProjection = True
        
        self.vertical_dim = None
        self.has_vertical_dim = False
        
        self.flipVerticalDim = False        
        self.verticalFactor = 0.0
        
        # these are for the particles which have missing values during the first time step(s) 
        # - eg. when the particle is starting later than other particles and so 
        # the lon, lat, vertical_dim array have NaN (missing values) at the beginning for the corresponding number of time steps
        self.stop_index = 1
        self.has_initial_missvals = False
        self.stop_indices = []
        self.first_timestep_missvals_indices = []
        
        # the variable holding the trajectory IDs
        # will have the attribute cf_role = "trajectory_id"
        self.nc_trajectory_id_var = None

        self.nc_point_variables = []

        #self.nc_dimensions = []
        self.timesteps = None 

        self.timeStepsCount = 0
        # will be possibly used in the future to figure out wether the timesteps are played forward or backward
        self.timeCompare = 0 



        from vtkmodules.vtkCommonCore import vtkDataArraySelection
        self.arrayselection = vtkDataArraySelection()
        self.arrayselection.AddObserver("ModifiedEvent", createModifiedCallback(self))


    def build_particlesID_array(self, numberOfParticles, numberOfTimesteps):
        vtk_array = vtk.vtkShortArray() 
        vtk_array.SetName("Particles ID")
        vtk_array.SetNumberOfComponents(1)
        vtk_array.SetNumberOfTuples(numberOfParticles*numberOfTimesteps)
        
        count = 0
        for i in range(0, numberOfTimesteps):
            for j in range(0, numberOfParticles):
                vtk_array.SetValue(count,i)
                count += 1

        return vtk_array
        

        
    
    def copy_first_valid_point(self):
        if self.has_initial_missvals: 
            # for each of the first time step NaNs, the first time step with a corresponding valid value is in self.stop_indices
            # we copy the first valid value over the missing values from the first time steps,
            # so that ParaView will alocate Pathlines when applying the Temporal Particles To Pathlines filter
            for k in range(0, len(self.first_timestep_missvals_indices)):
                index = self.first_timestep_missvals_indices[k]
                index_vert = self.stop_indices[k]

                for p in range(0, index_vert):
                    self.lat[p][index] = self.lat[index_vert][index]
                    self.lon[p][index] = self.lon[index_vert][index]
                    if self.has_vertical_dim == True: 
                        self.vertical_dim[p][index] = self.vertical_dim[index_vert][index]
            
                for m in self.nc_point_variables:
                    count = 0
                    for r in range(0, index_vert):
                        m.SetValue(index + count*self.numPoints, m.GetValue(index_vert*self.numPoints + index) ) 
                        count += 1
 



    def build_nc_array(self, nc_var):
        vtk_array = None

        if nc_var.dtype.type is np.float64:
            vtk_array = vtk.vtkDoubleArray()
        if nc_var.dtype.type is np.float32:
            vtk_array = vtk.vtkFloatArray()
        if nc_var.dtype.type is np.int8 or nc_var.dtype.type is np.int16:
            vtk_array = vtk.vtkShortArray()             
        if nc_var.dtype.type is np.int32 or nc_var.dtype.type is np.int64:
            vtk_array = vtk.vtkIntArray()

        #extract actual data from the data set         
        nc_var_data = nc_var[:]
        
        vtk_array.SetName(nc_var.name)
        vtk_array.SetNumberOfComponents(1)  
        vtk_array.SetNumberOfTuples(nc_var_data.size)         
        vtk_array.SetArray(nc_var_data.data, nc_var_data.size, True)
        vtk_array.array = nc_var_data.data
 
        return vtk_array


    def get_nc_data(self, requested_time=None):
  
        # get_nc_data will be called on each time step by RequestData - so make sure the data from the file is read in just once.
        if self.nc_dataset is not None:
            return self.nc_dataset

        if self.filename is None:
            # Note, exceptions are totally fine!
            raise RuntimeError("No filename specified")

        #read in the data
        self.nc_dataset = Dataset(self.filename, 'r+')

        # find and load the mandatory dimensions (time, lon, lat, trajectory_id and the generic variables of the dataset
        nc_info_var = [var for var in self.nc_dataset.variables.keys()] 
        nc_val_dim = self.nc_dataset.dimensions.values() 
        
        cf_role_match = ['cf_role']
        standard_match = ['standard_name']
        coordinates_match = ['coordinates']
        positive_match = ['positive']
        
        for i in nc_info_var: 
            nc_attribs = self.nc_dataset.variables[i].ncattrs()  

            cf_matching = [s for s in nc_attribs if any(xs in s for xs in cf_role_match)]
            
            # let's see which is the variable with the trajectory IDs 
            if len(cf_matching)>0 and cf_matching[0] == 'cf_role':
                if self.nc_dataset.variables[i].getncattr(cf_matching[0]) == 'trajectory_id':
                    print("The particles trajectories dimension: ", self.nc_dataset.variables[i].name) 
                    self.nc_trajectory_id_var = self.nc_dataset.variables[i][:]                    
                    
                    # let's identify the number of points from the corresponding particles trajectories dimension  
                    traj_dim_match = [self.nc_dataset.variables[i].name]
                    
                    for dimobj in nc_val_dim:
                        traj_matching = [s for s in [dimobj.name] if any(xs in s for xs in traj_dim_match)]  

                        if len(traj_matching)>0 and traj_matching[0] == traj_dim_match[0]:
                            #this is the trajectory dimension, it's size is the mumber of points (i.e. particles)
                            self.numPoints = dimobj.size
                           

            """
            A vertical coordinate will be identifiable by:          
                units of pressure; or       
                the presence of the positive attribute with a value of up or down (case insensitive).
                
             The Z  dimension will be identified by the "positive" attribute 
             If there is no Z dimension, we put a 0 value as the 3rd coordinate
            """
            pos_matching = [s for s in nc_attribs if any(xs in s for xs in positive_match)]  
            
            if len(pos_matching)>0 and pos_matching[0] == 'positive':       
                #this is the vertical coordinate dimension  
                self.vertical_dim = self.nc_dataset.variables[i][:]
                self.has_vertical_dim = True
                
                if self.nc_dataset.variables[i].getncattr(pos_matching[0]) == 'down':
                    # the vertical dimension is oriented "down", so its values will be made negative if necessary 
                    self.flipVerticalDim = True

                

            # let's see which variables have the "coordinates" attribute - this will be the actual variables 
            coords_matching = [s for s in nc_attribs if any(xs in s for xs in coordinates_match)]       
            if len(coords_matching)>0 and coords_matching[0] == 'coordinates':        

                # point data
                if self.nc_dataset.variables[i].ndim == 2: 
                    self.nc_point_variables.append(self.build_nc_array(self.nc_dataset.variables[i]))


            matching = [s for s in nc_attribs if any(xs in s for xs in standard_match)]

            if len(matching)>0 and matching[0] == 'standard_name':            
                if self.nc_dataset.variables[i].getncattr(matching[0]) == 'time':
                    self.time = self.nc_dataset.variables[i][:]
                elif self.nc_dataset.variables[i].getncattr(matching[0]) == 'latitude':
                    self.lat = self.nc_dataset.variables[i][:]
                    
                    """ 
                    We'll have to find out if any of the particles start with missing values
                    - for example when some particles start to move later, say only with the 5th time step, the values for time steps 1 to 4 will be filled wit missing values or NaNs
                    ParaView will NOT allocate Pathlines for them, when the Temporal Particles to Pathlines filter will be applied on the data.
                    So we'll have to find the indices of the first valid value for such particles, append them to self.stop_indices
                    then later copy the first valid value over those specific missing ones - this will be done with self.copy_first_valid_point()
                    This will result in ParaView allocation the correct number of Pathlines and Particles - and those specific Particles will 
                    appear directly on the first valid position in the visualization - adn will start to move accordingly, when the valid time steps are reached.
                    """
                    lat_missvals_indices = np.argwhere(ma.getmask(self.lat))
                    mask_lat = ma.getmask(self.lat) 
                    mask_lat_transpose = np.transpose(mask_lat)
                    
                    if len(lat_missvals_indices)> 0:
                        # the 0 element we already know it's a missval   --- Hmmm 
                        #self.first_timestep_missvals_indices.append(0) 
                                                
                        #This means there is some particle with a missing value at time step 0 
                        if lat_missvals_indices[0][0] == 0: 
                            
                            self.has_initial_missvals = True

                            for k in range(0, len(lat_missvals_indices)):
                                if lat_missvals_indices[k][0] == 0:
                                    self.first_timestep_missvals_indices.append(lat_missvals_indices[k][1])
                            
                            for c in range(0, len(self.first_timestep_missvals_indices) ) :
                                a = self.first_timestep_missvals_indices[c]

                                #look for index of first valid element - i.e. where there is no masked element (NaN or missing value)
                                v = np.argmax(mask_lat_transpose[a] == False)
                                self.stop_indices.append(v)
                    
                elif self.nc_dataset.variables[i].getncattr(matching[0]) == 'longitude':
                    self.lon = self.nc_dataset.variables[i][:]
                    
 

        # check if we have all four mandatory dimensions
        assert self.time is not None, "The time dimension is missing! The data will not be displayed correctly."
        assert self.lat is not None, "The latitude dimension is missing! The data will not be displayed correctly."
        assert self.lon is not None, "The longitude dimension is missing! The data will not be displayed correctly."        
        assert self.nc_trajectory_id_var is not None, "The particle trajectories variable is missing! The data will not be displayed correctly."


        self.timeStepsCount = len(self.time.data)

        # if we didn't find a vertical dimension while iterating above through nc_info_var, we'll set up an array of zeros so that all Z coords will be 0
        if self.vertical_dim is None: 
            self.vertical_dim = np.zeros((self.timeStepsCount,self.numPoints)) 
   
            
        """    
        As explained a few rows above, this is a fix for data sets that have *at the beginning* missing values:   
        Some particles may start not at the first time step, but at a later time step - so that they'd have missing values for some of the time steps *at the beginning*. 
        The "Temporal Pathlines" ParaView filter looks only at the *first* time step to allocate Particles Paths and will skip the points with missing values. 
        This is why, we will search the data set for the first valid particle position and values and will copy them over all the missing values from the corresponding time steps at the beginning.  
        """
        if self.has_initial_missvals:
            self.copy_first_valid_point()

        
        self.timesteps = None

        if len(self.time) > 0:
            self.timesteps = self.time

        for i in self.nc_point_variables:
            self.arrayselection.AddArray(i.GetName())
           

        return self.get_nc_data(requested_time)

    
    
    def get_timesteps(self):
        self.get_nc_data()
        return self.timesteps.tolist() if self.timesteps is not None else None

    
    
    def get_update_time(self, outInfo):
        executive = self.GetExecutive()
        timesteps = self.get_timesteps()
        if timesteps is None or len(timesteps) == 0:
            return None
        elif outInfo.Has(executive.UPDATE_TIME_STEP()) and len(timesteps) > 0:
            utime = outInfo.Get(executive.UPDATE_TIME_STEP())
            dtime = timesteps[0]
            for atime in timesteps:
                if atime > utime:
                    return dtime
                else:
                    dtime = atime
            return dtime
        else:
            assert(len(timesteps) > 0)
            return timesteps[0]

    def _get_array_selection(self):
        return self.arrayselection


    
    def buildPolyData(self, currentAnimationTime):    


        x = self.lon[currentAnimationTime]
        y = self.lat[currentAnimationTime]
        z = self.vertical_dim[currentAnimationTime]
        
        # the vertical dimension has the "positive" attribute "down", so it will be made negative if necessary:
        if self.flipVerticalDim:
            z = -abs(z)

        z_div = z/6370.0   #Vertical dimension scaled by the Earth radius in km
        
        if self.sphericalProjection:
            
            xr = np.radians(x) # LON
            yr = np.radians(y) # LAT
           
            # The unit sphere radius
            R = 1
            spherical_multip_coef = R + z_div*self.verticalFactor  

            xc = spherical_multip_coef*np.cos(yr)*np.cos(xr)  
            yc = spherical_multip_coef*np.cos(yr)*np.sin(xr)
            zc = spherical_multip_coef*np.sin(yr) 

            x = xc
            y = yc
            z = zc
        else: 
            z = z_div * self.verticalFactor

        if self.numPoints > 1 : 
            for j in range(0, self.numPoints):
                self.newPoints.InsertPoint(j, x[j], y[j], z[j])
        else:              
            self.newPoints.InsertPoint(self.polyIndex, x, y, z)  #they will aways be inserted at 0
    
        # Update the output polydata
        self.pdo.SetPoints(self.newPoints)        
        
        self.timeCompare = currentAnimationTime
     


    # GUI
    @smproperty.stringvector(name="FileName")
    @smdomain.filelist()
    @smhint.filechooser(extensions="nc", file_description="NetCDF particles trajectories files")
    def SetFileName(self, name):
        """Specify filename for the file to read."""
        if self.filename != name:
            self.filename = name
            self.nc_dataset = None
            self.timesteps = None
            self.Modified()

    @smproperty.doublevector(name="TimestepValues", information_only="1", si_class="vtkSITimeStepsProperty")
    def GetTimestepValues(self):
        return self.get_timesteps()

    #The VTK array selection API allows users to choose which arrays to
    #load. smproperty.dataarrayselection() exposes this in ParaView
    #This method *must* return a `vtkDataArraySelection` instance
    @smproperty.dataarrayselection(name="Arrays")
    def GetDataArraySelection(self):
        return self._get_array_selection()
        


    @smproperty.xml("""
       <IntVectorProperty name="SphericalProjection"
           number_of_elements="1"
           default_values="0"
           command="DisplaySphericalProjection">
           <BooleanDomain name="bool"/>
           <Documentation>Check this box for spherical projection.</Documentation>
       </IntVectorProperty>""")
    def DisplaySphericalProjection(self, x):
        if x == 1:
            self.sphericalProjection = True
        else:
            self.sphericalProjection = False
        self.Modified() 
      

    
    @smproperty.xml("""
       <DoubleVectorProperty  name="VerticalFactor"
           number_of_elements="1"
           default_values="1.0"
           command="SetVerticalFactor">
           <DoubleRangeDomain name="range" min="0.0" max="100000.0"/>
           <Documentation>This factor will be multiplied with the vertical dimension. Use the slider or type in a value  </Documentation>
       </DoubleVectorProperty >""")       
    #@smdomain.doublerange(min=0, max=100000.0)
    def SetVerticalFactor(self, x):        
        self.verticalFactor = x 
        self.Modified()


    # RequestInformation is called for the initial loading of the data
    def RequestInformation(self, request, inInfoVec, outInfoVec):

        executive = self.GetExecutive()
        outInfo = outInfoVec.GetInformationObject(0)
        outInfo.Remove(executive.TIME_STEPS())
        outInfo.Remove(executive.TIME_RANGE())

        #get data information and do initial data loading
        timesteps = self.get_timesteps()  

        if timesteps is not None:
            for t in timesteps:
                outInfo.Append(executive.TIME_STEPS(), t)
            outInfo.Append(executive.TIME_RANGE(), timesteps[0])
            outInfo.Append(executive.TIME_RANGE(), timesteps[-1])

            #allocate number of cells that will be added to the polyline                      
            self.pdo.Allocate(self.timeStepsCount*self.numPoints, 1)


        return 1



    #RequestData is called for the initial data loading and also whenever advancing the timesteps
    def RequestData(self, request, inInfoVec, outInfoVec):

        data_time = self.get_update_time(outInfoVec.GetInformationObject(0)) 

        #the array index of the currently requested time
        currentAnimationTime = np.where(self.time == data_time)
    
        self.get_nc_data(data_time)
        self.buildPolyData(currentAnimationTime[0][0]) 

        output = vtk.vtkPolyData.GetData(outInfoVec, 0)
        output.ShallowCopy(self.pdo)

        #populate the output with the point data from the current time step - i.e. extract the slice corresponding to the current time step
        for i in self.nc_point_variables:

            crt_array_slice = i.array[currentAnimationTime[0][0]][:]

            np_arr = np.array(crt_array_slice)
            
            vtk_crt_slice_array = None
            
            if np_arr.dtype.type is np.float64:
                vtk_crt_slice_array = vtk.vtkDoubleArray()
            if np_arr.dtype.type is np.float32:
                vtk_crt_slice_array = vtk.vtkFloatArray()
            if np_arr.dtype.type is np.int8 or np_arr.dtype.type is np.int16:
                vtk_crt_slice_array = vtk.vtkShortArray()             
            if np_arr.dtype.type is np.int32 or np_arr.dtype.type is np.int64:
                vtk_crt_slice_array = vtk.vtkIntArray()

            vtk_crt_slice_array.SetName(i.GetName())
            vtk_crt_slice_array.SetNumberOfComponents(1)
            vtk_crt_slice_array.SetNumberOfTuples(self.numPoints)  
            vtk_crt_slice_array.SetArray(np_arr, self.numPoints, True)
            vtk_crt_slice_array.array = np_arr
          
            output.GetPointData().AddArray(vtk_crt_slice_array)

        if data_time is not None:
            self.pdo.GetInformation().Set(self.pdo.DATA_TIME_STEP(), data_time)

        return 1
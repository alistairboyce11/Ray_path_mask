# Ray_path_mask
Use stat evt locs to calculate ray pierce points for approximate tomographic resolution masks

mk_ray_path_mask.py :

    def calc_mantle_fresnel_zones()
        Computes depth array of fresnel zone widths for P-waves.
    
    def cal_pierce_points()
        Takes AARM derived phase file input (Boyce et al., 2017; BSSA)
        Loops through all rays:
            where rays are not horizontal < tunring-depth*0.8:
                finds all pierce points for each ray at 50km depth intervals in mantle
                prints to file
            

    def haversine()
        Calculates distances between two points at depth in earth.
        Pass a point and a grid (vector of lats and lons points)
    
    def create_mask_grid()
        Takes Ray PP file generated above.
        Computes a grid equal to tomographic plotting grid
        for each depth finds all ray PPs:
            for each PP:
                Finds distance from PP to each point in grid.
                Marks as zero if less than fzhw*fzmf+fzgs
                (fzmf = 'Fresnel zone multiplication factor' Used to broaden the footprint of each fresnel zone based mask.
                (fzgs = 'Fresnel zone grid size' A static value added to region "captured" by ray. Tomography is gridded in constant velocity cells of minimum size ~fzgs)
                sums all PP grid at each depth to finalise each mask depth interval.
            
parallel_example.py :

    Example script to practice code parallelisation.

par_mk_ray_mask.py : 

    This is parallelized code for computing mask grids from the PP files.
    Not written in nested functions as cannot be pickled for passage to/from cores.
    Must change params at top of script
    As above for create_mask_grid(), but:
        for each PP concurrent.futures.ProcessPoolExecutor takes care of computing the distances to the grid at each depth interval.
        Massively speed up as avoids a for-loop greater than length of ray pp infile.
        
    

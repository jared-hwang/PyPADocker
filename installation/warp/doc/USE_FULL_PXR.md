When PICSAR is compiled using "full_pxr" mode (See PICSAR documentation), WARP can take advantage
of additional features from picsar. This notably includes the Hybrid PSATD solver capability for solving Maxwell's equations as well as absorbing boundary conditions.
This mode allows a substantial memory saving due to overall less data redundancy

Hybrid PSATD allows to solve Maxwell's equations using distributed memory FFT across multiple mpi subdomains in order to reduce the memory footprint and the computational time of the simulation. At present, PICSAR supports
P3DFFT and FFTW libraries for performing distributed FFTs. The decomposition technique used groups the mpi tasks among different "mpi groups", each mpi group still forming a cartesian subdomain.
Then distributed ffts are performed across each mpi group and guard cells exchanged between the mpi groups. This results in better memory footprint (tanks to reduced data redundancy) as well as better performance.

Note that when using full_pxr = True, Perflectly Matched Layer (PML) computations for simulating absorbing boundary conditions are done in PICSAR exclusively. In PICSAR, PMLs are inside the simulation domain you set, (the first and last np_pml cells are set to be PML cells)
You might take that into consideration when designing you simulation box.


When using this mode, you need to turn the `full_pxr` flag in your EM3DPXR class instanciation to `True`.
Then you need to specify additional parameters:


em=EM3DPXR{.
	   .
	   .
           'full_pxr': False,
           'fftw_hybrid':False,
           'fftw_with_mpi':False,
           'fftw_mpi_transpose':False,
           'p3dfft_flag':False,
           'p3dfft_stride':False,
	   'shift_x_pml_pxr':4,
           'shift_y_pml_pxr':4,
           'shift_z_pml_pxr':4,
           'nb_group_x':0,
           'nb_group_y':0,
           'nb_group_z':0,
           'nyg_group':0,
           'nzg_group':0,
	   'absorbing_bcs_x':False,
           'absorbing_bcs_y':False,
           'absorbing_bcs_z':False,
           'nx_pml':8,
           'ny_pml':8,
           'nz_pml':8,
           'g_spectral':False,
           }

when `full_pxr` is set to `True`, this will allow matrix initializations for PSATD and PMLs compuations to be done in PICSAR (initialization on the WARP layer is not performed).
To use the hybrid PSATD solver, you need to set `fftw_with_mpi` to `True` and `fftw_hybrid=True`. If `fftw_hyrid = False`, then PICSAR will solve Maxwell's equations 
using a global domain FFT. Note that this mode is buggy and has not been tested extensively. As a consequence, use only `fftw_with_mpi = True`, and `fftw_hybrid = True` for now in order to avoid inconvenient bugs.

Regarding the choice of the FFT library, you ca choose between FFTW and P3DFFT. To choose P3DFFT set `p3dfft_flag = True`, otherwise warp will use FFTW_MPI by default.

Whether you chose FFTW or P3DFFT, you might want to set `fftw_mpi_transpose = True` or `p3dfft_stride = True` respectively. This performs faster FFTs since it gets rid of one data transposition during the distributed FFT computation.
Note that when using `p3dfft_stride = True`, the P3DFFT library needs to be compiled using `--enable-stride1`. In addition, when `p3dfft_stride=False`, you need to recompile P3DFFT without `--enable-stride1`. However, we recommend to always use `--enable-stride1` since this results in better performance.


The main difference between P3DFFT and FFTW libraries is that FFTW allows for 1D slab domain decomposition (along z-axis) whereas P3DFFT allows for 2D pencil decomposition. Note that P3DFFT CANNOT be used when doing 2D simulations since this library only performs 3D FFTs.
Plus, `fftw_mpi_transpose` must be set to `False` when performing 2D simulations.


The next parameter to be set are `nb_group_x,y,z`: 

`nb_group_x` is the number of groups along the x direction.
`nb_group_y` is the number of groups along the y direction.
`nb_group_z` is the number of groups along the z direction.

At present, you can't choose `nb_group_x` as there is currently no support for FFT libraries allowing 3D decomposition. If you are using FFTW (1D slab decomposition), nb_group_y will be set automatically to nprocy.  

You will also need to set `nyg_group`(number of guardcells along y for the groups)
and `nzg_group`(number of guardcells along z  for the groups). When using the hybrid mode, you can set `nzg_group`, `nyg_group` (P3DFFT only) high enough to get very small truncation errors and reduce `nyguards`, `nzguards` to the minimum value to avoid data redundancy. 
`nxg_group` is set to `nxguards` automatically since task gathering along the x-axis is currently not supported.


If you are doing simulations using periodic boundary counditions for the EM-fields you might want to set `g_spectral = False` as this will reduce the memory footprint of the EM solver.

When using full_pxr mode with pmls you need to set absorbing_bcs_x/y/z flags to True. Unlike warp, picsar can handle different EM boundary conditions along y and z (but same BC for z_min and z_max).
You would also want force w3d boundaries on fields to False.

In general you don't need to change the values of `nx_pml`, `ny_pml` ... since the default value already brings good absorbtion rate.
The PML in picsar in by default appended to the first  and the last ni_pmls OF the simulation domain. Field values inside domain guardcells are forced to 0 to emulate perfectly reflective media.
When using the local solver with full pxr this can be changed slightly in order to put the pml cells inside the domain guardcells though the shift_x_pml/shift_y_pml/shift_z_pml parameters.
These parameters represent the number of domain guardcells in which the fields are forced to 0. So if nx_pml + shift_x_pml = nxguards for example then the whole pml region is held by the guardcell.
For hybird solver, this has not been implemented yet, and shift parameters don't have any impact









##############################################################################
# IA32_Intel/EM64T_Intel:	Intel Pentium with MS Windows using Intel compiler 11.
##############################################################################
FLAGS_OPT        = -O    
FLAGS_OMP        =  -fopenmp
H2Ocean_EXE      = h2ocean
EXTO             = o 
RM               = rm
FOR              = gfortran
FLINKER          = gfortran
INCS             =
LIBS             =
MAKE             = make
.SUFFIXES: .f .f90 .F90 .o .obj


OBJECTS = \
Mod_Precision.$(EXTO) \
Mod_Sparskit.$(EXTO) \
Mod_Sparse_Matrix.$(EXTO) \
Mod_bivar.$(EXTO) \
Mod_Parameter.$(EXTO) \
Mod_IO.$(EXTO) \
Mod_Derived_DataType.$(EXTO) \
Mod_All_Variable.$(EXTO) \
Mod_TriMesh.$(EXTO) \
Mod_TriMesh_Interp.$(EXTO) \
Mod_FEM.$(EXTO) \
Mod_H2Ocean.$(EXTO) \
Mod_Control_File.$(EXTO) \
Mod_Dirichlet_Boundary.$(EXTO) \
Subs_Ctrl_File.$(EXTO) \
H2Ocean_Version_Info.$(EXTO) \
H2Ocean_Read_Mesh_File.$(EXTO) \
Subs_H2Ocean_Time_Stamp.$(EXTO) \
Subs_H2Ocean_Initialize.$(EXTO) \
Subs_H2Ocean_Output.$(EXTO) \
Subs_H2Ocean_alpha_h.$(EXTO) \
Subs_H2Ocean_Compute_dt.$(EXTO) \
Subs_H2Ocean_Nonhydro.$(EXTO) \
Subs_H2Ocean_Computation.$(EXTO) \
Subs_H2Ocean_Finalize.$(EXTO) \
H2Ocean_Main.$(EXTO)


ser:
	$(MAKE) FFLAGS="$(FLAGS_OPT)" \
	        FFLAGS90="$(FLAGS_OPT)" $(H2Ocean_EXE)	
omp:
	$(MAKE) FFLAGS="$(FLAGS_OPT)  $(FLAGS_OMP)" \
	        FFLAGS90="$(FLAGS_OPT)  $(FLAGS_OMP)" $(H2Ocean_EXE)

$(H2Ocean_EXE): $(OBJECTS)  
	-$(FLINKER) $(FFLAGS) -o  $(H2Ocean_EXE) $(OBJECTS) $(INCS) $(LIBS)

.PHONEY: clean

.f90.o:
	$(FOR) $< -c $(FFLAGS90) $(INCS)  
.f90.obj:
	$(FOR) $< -c $(FFLAGS90)  $(INCS)   
.F90.o:
	$(FOR) $< -c $(FFLAGS90)  $(INCS) 
.F90.obj:
	$(FOR) $< -c $(FFLAGS90)  $(INCS) 
.f.o:
	$(FOR) $< -c $(FFLAGS)   $(INCS)
.f.obj:
	$(FOR) $< -c $(FFLAGS)   $(INCS)
#---------------------------------------------------------------------------- 	
clean:	
	-@$(RM) *.$(EXTO); 
	-@$(RM) *.mod; 
	-@$(RM) $(H2Ocean_EXE);

install:
	h2ocean= `pwd`/h2ocean
	export h2ocean

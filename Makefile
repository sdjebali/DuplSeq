#######################
#    OCAML Programs   #
#######################

OCAMLC = ocamlc
OCAMLOPT = ocamlopt
OCAMLDEP = ocamldep



#######################
#     OCAML Flags     #
#######################

INCLUDES =					#all relevant -I options here
OCAMLFLAGS = $(INCLUDES)	#add other options for ocamlc here
OCAMLOPTFLAGS = $(INCLUDES)	#add other options for ocamlopt here



#######################
#    Sources files    #
#######################

MLFILES = common.ml config.ml feature.ml input.ml output.ml main.ml
CMXAFILES = $(CMAFILES:%.cma=%.cmxa) 
CMOFILES = $(MLFILES:%.ml=%.cmo) 
CMXFILES = $(MLFILES:%.ml=%.cmx)
CMIFILES = $(MLFILES:%.ml=%.cmi) $(MLIFILES:%.mli=%.cmi)
OBJFILES = $(CMIFILES:%.cmi=%.o)
BINFILE = duplseq_gal



#######################
#        Rules        #
#######################

duplseq_gal: $(CMXFILES)
	@echo "LNKOPT $(BINFILE)"
	@$(OCAMLOPT) $(CMXAFILES) $(CMXFILES) $(OCAMLOPTFLAGS) -o $(BINFILE)


# Common rules
.SUFFIXES: .ml .cmo .cmi .cmx

.ml.cmo:
	@echo "OCAMLC $<"
	@$(OCAMLC) $(OCAMLFLAGS) -c $<

.mli.cmi:
	@echo "OCAMLC $<"
	@$(OCAMLC) $(OCAMLFLAGS) -c $<

.ml.cmx:
	@echo "OCAMLOPT $<"
	@$(OCAMLOPT) $(OCAMLOPTFLAGS) -c $<

# Clean up
clean:
	@echo "Cleaning .cmo .cmx .o and .cmi"
	@rm -f $(CMOFILES) $(CMXFILES) $(OBJFILES) 

# Dependencies
depend:
	@echo "Calculating dependencies"
	@$(OCAMLDEP) $(INCLUDES) $(MLFILES) $(MLIFILES) > $(DEPENDFILE)

# include $(DEPENDFILE)



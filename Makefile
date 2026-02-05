#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*    This file is part of the code gasa01, which contains code              *
#*    for Project A01 "Global Methods for Stationary Gastransport"           *
#*    of the SFB#TRR 154 "Mathematical Modelling, Simulation and             *
#*    Optimization on the Example of Gas Networks".                          *
#*                                                                           *
#*    Copyright (C) 2015-     Discrete Optimization Group, TU Darmstadt      *
#*                                                                           *
#*    Based on SCIP  --- Solving Constraint Integer Programs                 *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#@file    Makefile
#@brief   Makefile for gasa01
#@author  Marc Pfetsch


#-----------------------------------------------------------------------------
# paths
#-----------------------------------------------------------------------------

SCIPDIR         =       $(SCIP_PATH)
SCIPREALPATH	=	$(realpath $(SCIPDIR))

#-----------------------------------------------------------------------------
# include default project Makefile from SCIP
#-----------------------------------------------------------------------------

include $(SCIPDIR)/make/make.project

VERSION		=	2.0
TEST		=	basic
DFLAGS		= 	-MM

#-----------------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------------

MAINNAME	=	gasa01
MAINOBJ		=	cmain.o

MAINSRC		=	$(addprefix $(SRCDIR)/,$(MAINOBJ:.o=.c))
MAINDEP		=	$(SRCDIR)/depend.cmain.$(OPT)

MAIN		=	$(MAINNAME).$(BASE).$(LPS)$(EXEEXTENSION)
MAINFILE	=	$(BINDIR)/$(MAIN)
MAINSHORTLINK	=	$(BINDIR)/$(MAINNAME)
MAINOBJFILES	=	$(addprefix $(OBJDIR)/,$(MAINOBJ))

GASOBJ		=	probdata_gas.o generateModel_gas.o unitConversion_gas.o characteristics_gas.o ode_functions_gas.o \
			cons_pipe_ode.o LSFiles_gas.o disp_gradient.o comparisons_gas.o heur_gas.o presol_flowobbt.o heur_ADM.o heur_lpround.o
GASSRC		=	$(addprefix $(SRCDIR)/,$(GASOBJ:.o=.c))
GASOBJFILES	=	$(addprefix $(OBJDIR)/,$(GASOBJ))

#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

ifeq ($(VERBOSE),false)
.SILENT:	$(MAINFILE) $(MAINOBJFILES) $(GASOBJFILES) $(MAINSHORTLINK) $(BINDIR)/test_ode_approx test_ode_approx $(OBJDIR)/test_ode_approx.o
endif

.PHONY: all
all:            $(SCIPDIR) $(MAINFILE) $(MAINSHORTLINK)

.PHONY: lint
lint:		$(MAINSRC) $(GASSRC)
		-rm -f lint.out
		$(SHELL) -ec 'for i in $^; \
			do \
			echo $$i; \
			$(LINT) lint/main-gcc.lnt +os\(lint.out\) -u -zero \
			$(USRFLAGS) $(FLAGS) -I/usr/include -UNDEBUG -USCIP_WITH_READLINE -USCIP_ROUNDING_FE -D_BSD_SOURCE $$i; \
			done'

.PHONY: pclint
pclint:		$(MAINSRC) $(GASSRC)
		-rm -f pclint.out
ifeq ($(FILES),)
		@$(SHELL) -ec 'echo "-> running pclint ..."; \
			for i in $^; \
			do \
				echo $$i; \
				$(PCLINT) -I$(SCIPREALPATH) -I$(SCIPREALPATH)/pclint main-gcc.lnt +os\(pclint.out\) -b -u -zero \
				$(USRFLAGS) $(FLAGS) $(SDPIINC) -uNDEBUG -uSCIP_WITH_READLINE -uSCIP_ROUNDING_FE -D_BSD_SOURCE $$i; \
			done'
else
		@$(SHELL) -ec  'echo "-> running pclint on specified files ..."; \
			for i in $(FILES); \
			do \
				echo $$i; \
				$(PCLINT) -I$(SCIPREALPATH) -I$(SCIPREALPATH)/pclint main-gcc.lnt +os\(pclint.out\) -b -u -zero \
				$(USRFLAGS) $(FLAGS) $(SDPIINC) -uNDEBUG -uSCIP_WITH_READLINE -uSCIP_ROUNDING_FE -D_BSD_SOURCE $$i; \
			done'
endif


.PHONY: scip
scip:
		@$(MAKE) -C $(SCIPDIR) libs $^

.PHONY: doc
doc:
		cd doc; $(DOXY) $(MAINNAME).dxy

$(MAINSHORTLINK):	$(MAINFILE)
		@rm -f $@
		cd $(dir $@) && ln -s $(notdir $(MAINFILE)) $(notdir $@)

$(OBJDIR):
		@-mkdir -p $(OBJDIR)

$(BINDIR):
		@-mkdir -p $(BINDIR)

.PHONY: clean
clean:		$(OBJDIR)
ifneq ($(OBJDIR),)
		@-(rm -f $(OBJDIR)/*.o && rmdir $(OBJDIR));
		@echo "-> remove main objective files"
endif
		@-rm -f $(MAINFILE) $(MAINLINK) $(MAINSHORTLINK)
		@echo "-> remove binary"

.PHONY: test
test:           $(MAINFILE)
		cd check; \
		$(SHELL) ./check.sh $(TEST) $(MAINFILE) $(SETTINGS) $(notdir $(MAINFILE)) $(TIME) $(NODES) $(MEM) $(DISPFREQ) $(VERSION) $(LPS);

.PHONY: testcluster
testcluster:	$(MAINFILE)
		cd check; \
		$(SHELL) ./check_cluster.sh $(TEST) $(MAINFILE) $(SETTINGS) $(notdir $(MAINFILE)) $(TIME) $(NODES) $(MEM) $(DISPFREQ) $(VERSION) \
			$(LPS) $(QUEUE) $(CLIENTTMPDIR) $(OPT);

.PHONY: testlocal
testlocal: $(MAINFILE)
		cd check; \
		$(SHELL) ./check_local.sh $(TEST) $(MAINFILE) $(SETTINGS) $(TIME) $(NODES) $(MEM) $(DISPFREQ) $(VERSION) $(CORES) $(LPS) ;

.PHONY: tags
tags:
		rm -f TAGS; etags src/*.c src/*.h $(SCIPDIR)/src/*/*.c $(SCIPDIR)/src/*/*.h --output=TAGS; sed 's!\#undef .*!!g' TAGS > tags; mv tags TAGS


.PHONY: depend
depend:		$(SCIPDIR)
		$(SHELL) -ec '$(DCC) $(FLAGS) $(DFLAGS) $(MAINSRC) $(GASSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/\([0-9A-Za-z\_]*\).c|$$\(OBJDIR\)/\2.o: $(SRCDIR)/\2.c|g'\'' \
		>$(MAINDEP)'

-include	$(MAINDEP)


# main target
$(MAINFILE):	$(BINDIR) $(OBJDIR) $(SCIPLIBFILE) $(LPILIBFILE) $(NLPILIBFILE) $(MAINOBJFILES) $(GASOBJFILES)
		@echo "-> linking $@"
ifdef LINKCXXSCIPALL
		$(LINKCXX) $(GASOBJFILES) $(MAINOBJFILES) $(LINKCXXSCIPALL) $(LINKCXX_o)$@
else
		$(LINKCXX) $(GASOBJFILES) $(MAINOBJFILES) \
		$(LINKCXX_L)$(SCIPDIR)/lib $(LINKCXX_l)$(SCIPLIB)$(LINKLIBSUFFIX) \
                $(LINKCXX_l)$(OBJSCIPLIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(LPILIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(NLPILIB)$(LINKLIBSUFFIX) \
                $(OFLAGS) $(LPSLDFLAGS) \
		$(LDFLAGS) $(LINKCXX_o)$@
endif

# test target
test_ode_approx: $(BINDIR)/test_ode_approx
$(BINDIR)/test_ode_approx: $(BINDIR) $(OBJDIR) $(SCIPLIBFILE) $(LPILIBFILE) $(NLPILIBFILE) $(GASOBJFILES) $(OBJDIR)/test_ode_approx.o
		@echo "-> linking $@"
ifdef LINKCXXSCIPALL
		$(LINKCXX) $(GASOBJFILES) $(OBJDIR)/test_ode_approx.o $(LINKCXXSCIPALL) $(LINKCXX_o)$@
else
		$(LINKCXX) $(GASOBJFILES) $(OBJDIR)/test_ode_approx.o \
		$(LINKCXX_L)$(SCIPDIR)/lib $(LINKCXX_l)$(SCIPLIB)$(LINKLIBSUFFIX) \
                $(LINKCXX_l)$(OBJSCIPLIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(LPILIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(NLPILIB)$(LINKLIBSUFFIX) \
                $(OFLAGS) $(LPSLDFLAGS) \
		$(LDFLAGS) $(LINKCXX_o)$@
endif

$(OBJDIR)/%.o:	$(SRCDIR)/%.c
		@echo "-> compiling $@"
		$(CC) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CFLAGS) -c $< $(CC_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@echo "-> compiling $@"
		$(CXX) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) -c $< $(CXX_o)$@

#---- EOF --------------------------------------------------------------------

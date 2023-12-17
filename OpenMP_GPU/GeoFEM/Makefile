subdirs=Omp_U_f_a Omp_U_f_c Omp_U_s_a Omp_U_s_c

.PHONY: all $(subdirs)

all: $(subdirs)

$(subdirs):
	$(MAKE) -C $@ $(MAKECMDGOALS)

clean: $(subdirs)


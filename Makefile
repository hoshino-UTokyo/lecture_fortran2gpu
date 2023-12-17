subdirs=DoConcurrent OpenACC OpenMP_CPU OpenMP_GPU OpenACC+CUDA

.PHONY: all $(subdirs)

all: $(subdirs)

$(subdirs):
	$(MAKE) -C $@ $(MAKECMDGOALS)

clean: $(subdirs)


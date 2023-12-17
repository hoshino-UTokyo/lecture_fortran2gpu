subdirs=CUDA_U_f_a

.PHONY: all $(subdirs)

all: $(subdirs)

$(subdirs):
	$(MAKE) -C $@ $(MAKECMDGOALS)

clean: $(subdirs)


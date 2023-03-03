

FC=gfortran

DASSL_FLAGS = -std=legacy

linpack.o: linpack/linpack.f
	$(FC) -c $<

.phony: dassl

dassl: ddaini.o ddajac.o ddanrm.o ddaslv.o ddassl.o \
	   ddastp.o ddatrp.o ddawts.o s88fmt.o scherr.o xerrwv.o


ddaini.o: dassl/ddaini.f
	$(FC) $(DASSL_FLAGS) -c $<

ddajac.o: dassl/ddajac.f
	$(FC) $(DASSL_FLAGS) -c $<

ddanrm.o: dassl/ddanrm.f
	$(FC) $(DASSL_FLAGS) -c $<

ddaslv.o: dassl/ddaslv.f
	$(FC) $(DASSL_FLAGS) -c $<

ddassl.o: dassl/ddassl.f
	$(FC) $(DASSL_FLAGS) -c $<

ddastp.o: dassl/ddastp.f
	$(FC) $(DASSL_FLAGS) -c $<

ddatrp.o: dassl/ddatrp.f
	$(FC) $(DASSL_FLAGS) -c $<

ddawts.o: dassl/ddawts.f
	$(FC) $(DASSL_FLAGS) -c $<

s88fmt.o: dassl/s88fmt.f
	$(FC) $(DASSL_FLAGS) -c $<

scherr.o: dassl/scherr.f
	$(FC) $(DASSL_FLAGS) -c $<

xerrwv.o: dassl/xerrwv.f
	$(FC) $(DASSL_FLAGS) -c $<

.phony: clean

clean:
	rm -rf *.o
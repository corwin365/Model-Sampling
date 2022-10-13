function  InstrumentTypes = oif_list();

%list of the output directories used by all the Observation Import Functions currently implemented.
%this is only used when settings up paths using setup_instructions.m
%so it's not critical to keep it up to date (but it is helpful!)

InstrumentTypes = {'AIRS','AIRS_1D','AIRS_3D','AIRS_qbo','COSMIC','COSMICpetr','gisinger','HIRDLS','PETERP','SABER'}

return

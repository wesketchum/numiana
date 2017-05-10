# basic source code CMakeLists.txt

art_make( LIBRARY_NAME numiana_ChargedTrackMultiplicity
	  	 MODULE_LIBRARIES
			nusimdata_SimulationBase
			lardataobj_MCBase
			art_Utilities canvas_Utilities
			art_Persistency_Common canvas_Persistency_Common
			art_Persistency_Provenance canvas_Persistency_Provenance
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_IO_SOURCES}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${FHICLCPP}
                        ${CETLIB}
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${ROOT_BASIC_LIB_LIST}
		)

install_headers()
install_source()
install_fhicl()

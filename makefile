release:
	cd zoning; rm -f zoning.exe; make; cd ..
	mkdir -p x64/release
	mv zoning/zoning.exe x64/release/zoning.exe
	
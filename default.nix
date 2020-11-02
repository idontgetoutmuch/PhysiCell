{ stdenv
, fetchurl
, llvmPackages
, overrideCC
, gcc }:

let

openmp = llvmPackages.openmp;

in (overrideCC stdenv gcc).mkDerivation rec {
  pname = "PhysiCell";
  version = "1.7.1";

  buildInputs = [ openmp ];
  nativeBuildInputs = [ openmp ];

  src = fetchurl {
    url = "https://github.com/MathCancer/${pname}/archive/${version}.tar.gz";
    sha256 = "000pkw8maldb25b8g9cgr5qi4l905073pqrr3nfazc9n990260fi";
  };

  PHYSICELL_CPP="${gcc}/bin/g++";

  project = "biorobots";

  buildPhase = ''
    make $project-sample
    make
  '';

  installPhase = ''
    mkdir -p "$out"/bin
    cp $project $out/bin/$project
  '';

  meta = with stdenv.lib; {
    description = "PhysiCell: an Open Source Physics-Based Cell Simulator for 3-D Multicellular Systems.";
    homepage    = http://physicell.org;
    platforms   = platforms.all;
    maintainers = with maintainers; [ idontgetoutmuch ];
    license     = licenses.bsd3;
  };
}

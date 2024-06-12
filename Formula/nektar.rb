class Nektar < Formula
  desc "High-performance spectral/hp element framework"
  homepage "https://www.nektar.info/"
  url "https://gitlab.nektar.info/nektar/nektar/-/archive/v5.6.0/nektar-v5.6.0.tar.bz2"
  sha256 "215f16b7d62149f49206d71302b36266471a2bd7e45f328cad8a47f16e425c73"

  depends_on "arpack"
  depends_on "boost"
  depends_on "boost-python3"
  depends_on "cmake"
  depends_on "fftw"
  depends_on "hdf5-mpi"
  depends_on "numpy"
  depends_on "open-mpi"
  depends_on "opencascade"
  depends_on "python-setuptools"
  depends_on "python@3.12"
  depends_on "scotch"
  depends_on "tinyxml"
  depends_on "zlib"

  # Patch for HDF 1.12
  patch :DATA

  def install
    args = std_cmake_args + ["-DNEKTAR_BUILD_TESTS=OFF",
                             "-DNEKTAR_BUILD_UNIT_TESTS=OFF",
                             "-DZLIB_ROOT=#{Formula["zlib"].opt_prefix}"]

    args << "-DNEKTAR_BUILD_DEMOS=OFF"
    args << "-DNEKTAR_BUILD_MESHGEN=ON"
    args << "-DNEKTAR_BUILD_PYTHON=ON"
    args << "-DNEKTAR_SOLVER_DIFFUSION=OFF"
    args << "-DNEKTAR_SOLVER_DUMMY=OFF"
    args << "-DNEKTAR_SOLVER_ELASTICITY=OFF"
    args << "-DNEKTAR_SOLVER_MMF=OFF"
    args << "-DNEKTAR_USE_MPI=ON"
    args << "-DNEKTAR_USE_ARPACK=ON"
    args << "-DNEKTAR_USE_FFTW=ON"
    args << "-DNEKTAR_USE_SCOTCH=ON"
    args << "-DNEKTAR_USE_ARPACK=ON"
    args << "-DNEKTAR_USE_HDF5=ON"
    args << "-DPYTHON_EXECUTABLE=#{Formula["python@3.12"].opt_bin}/python3.12"

    mkdir "build" do
      rm "../cmake/FindHDF5.cmake"
      system "cmake", "..", *args
      system "make", "install"

      # Also need to install NekPy bindings
      chdir "python" do
        python = Formula["python@3.12"]
        system python.bin/"python3.12", *Language::Python.setup_install_args(prefix)

        site_packages = Language::Python.site_packages("python3.12")
        pth_contents = "import site; site.addsitedir('#{libexec/site_packages}')\n"
        (prefix/site_packages/"homebrew-nektar.pth").write pth_contents
      end
    end
  end

  test do
    (testpath/"helm.cpp").write <<~EOS
      #include <iostream>
      #include <MultiRegions/ContField.h>
      #include <SpatialDomains/MeshGraph.h>

      using namespace Nektar;
      using namespace std;

      int main(int argc, char *argv[])
      {
          LibUtilities::SessionReaderSharedPtr vSession
              = LibUtilities::SessionReader::CreateInstance(argc, argv);
          SpatialDomains::MeshGraphSharedPtr graph =
              SpatialDomains::MeshGraph::Read(vSession);
          MultiRegions::ContFieldSharedPtr fld =
              MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(
                  vSession, graph, vSession->GetVariable(0));

          // Set up forcing
          const int nq = fld->GetNpoints(), nm = fld->GetNcoeffs();
          Array<OneD, NekDouble> xc0(nq), xc1(nq), xc2(nq);
          fld->GetCoords(xc0, xc1, xc2);
          vSession->GetFunction("Forcing", 0)->Evaluate(xc0, xc1, xc2, fld->UpdatePhys());

          // Solve Helmholtz
          StdRegions::ConstFactorMap factors;
          factors[StdRegions::eFactorLambda] = vSession->GetParameter("Lambda");
          Vmath::Zero(nm, fld->UpdateCoeffs(), 1);
          fld->HelmSolve(fld->GetPhys(), fld->UpdateCoeffs(), factors);

          // Evaluate L^2 error
          fld->BwdTrans(fld->GetCoeffs(), fld->UpdatePhys());
          Array<OneD, NekDouble> sol(nq);
          vSession->GetFunction("ExactSolution", 0)->Evaluate(xc0, xc1, xc2, sol);
          cout << fld->L2(fld->GetPhys(), sol) << endl;
          return 0;
      }
    EOS

    (testpath/"CMakeLists.txt").write <<~EOS
      set(CMAKE_CXX_STANDARD 17)
      find_package(Nektar++ REQUIRED NO_MODULE NO_DEFAULT_PATH NO_CMAKE_BUILDS_PATH NO_CMAKE_PACKAGE_REGISTRY)
      include_directories(${NEKTAR++_INCLUDE_DIRS} ${NEKTAR++_TP_INCLUDE_DIRS})
      add_definitions(${NEKTAR++_DEFINITIONS})
      link_directories(${NEKTAR++_LIBRARY_DIRS} ${NEKTAR++_TP_LIBRARY_DIRS})
      add_executable(helm helm.cpp)
      target_link_libraries(helm ${NEKTAR++_LIBRARIES} ${NEKTAR++_TP_LIBRARIES})
    EOS

    (testpath/"project.py").write <<~EOS
      from NekPy.LibUtilities import PointsKey, PointsType, BasisKey, BasisType
      from NekPy.StdRegions import StdQuadExp
      import numpy as np
      ptsKey   = PointsKey(9, PointsType.GaussLobattoLegendre)
      basisKey = BasisKey(BasisType.Modified_A, 8, ptsKey)
      quadExp  = StdQuadExp(basisKey, basisKey)
      x, y     = quadExp.GetCoords()
      fx       = np.sin(x) * np.cos(y)
      proj     = quadExp.FwdTrans(fx)
      assert np.allclose(fx, quadExp.BwdTrans(proj))
    EOS

    (testpath/"input.xml").write <<~EOS
      <NEKTAR>
          <GEOMETRY DIM="2" SPACE="2">
              <VERTEX COMPRESSED="B64Z-LittleEndian" BITSIZE="64">eJxjYEAGD/aj0gwMjAzYAEKeCas8AjDj0AcDLNjNt4exWLG7Dy7PhlUfwh52HObCAAd2/XB1AAerEJkA</VERTEX>
              <EDGE COMPRESSED="B64Z-LittleEndian" BITSIZE="64">eJxjYEAFjDhoJhw0Mw4aBljQ1MP4rDj4bGh8mHnsaO6BqeNA48PUcaLxYfZzoYnD9HOj8WHuAgBAqACe</EDGE>
              <ELEMENT>
                  <Q COMPRESSED="B64Z-LittleEndian" BITSIZE="64">eJxjYEAFjFCaCUoz4xBngdKsUJoNTZ4dSnNAaU40c5jRxLmgNDea+QAUwABZ</Q>
              </ELEMENT>
              <COMPOSITE>
                  <C ID="0"> Q[0-3] </C>
                  <C ID="1"> E[0,3,5,6,7,8,10,11] </C>
              </COMPOSITE>
              <DOMAIN> C[0] </DOMAIN>
          </GEOMETRY>
          <EXPANSIONS>
              <E COMPOSITE="C[0]" NUMMODES="9" TYPE="MODIFIED" FIELDS="u" />
          </EXPANSIONS>
          <CONDITIONS>
              <PARAMETERS> <P> Lambda = 1 </P> </PARAMETERS>
              <VARIABLES> <V ID="0"> u </V> </VARIABLES>
              <BOUNDARYREGIONS> <B ID="0"> C[1] </B> </BOUNDARYREGIONS>
              <BOUNDARYCONDITIONS>
                  <REGION REF="0">
                      <D VAR="u" VALUE="sin(PI*x)*sin(PI*y)" />
                  </REGION>
              </BOUNDARYCONDITIONS>
              <FUNCTION NAME="Forcing">
                  <E VAR="u" VALUE="-(Lambda + 2*PI*PI)*sin(PI*x)*sin(PI*y)" />
              </FUNCTION>
              <FUNCTION NAME="ExactSolution">
                  <E VAR="u" VALUE="sin(PI*x)*sin(PI*y)" />
              </FUNCTION>
          </CONDITIONS>
      </NEKTAR>
    EOS
    system "cmake", "-DNektar++_DIR=#{lib}/nektar++/cmake/", "."
    system "make"
    assert (`./helm input.xml`.to_f < 1.0e-8)
    python = Formula["python@3.12"]
    system python.bin/"python3.12", "project.py"
  end
end

__END__
diff --git a/cmake/NektarCommon.cmake b/cmake/NektarCommon.cmake
index 01274b4fe..e14263b46 100644
--- a/cmake/NektarCommon.cmake
+++ b/cmake/NektarCommon.cmake
@@ -291,8 +291,14 @@ MACRO(ADD_NEKPY_LIBRARY name)
     TARGET_LINK_LIBRARIES(_${name}
         ${Boost_SYSTEM_LIBRARY}
         ${BOOST_PYTHON_LIB}
-        ${BOOST_NUMPY_LIB}
-        ${PYTHON_LIBRARIES})
+        ${BOOST_NUMPY_LIB})
+
+    IF(APPLE)
+        SET_TARGET_PROPERTIES(_${name} PROPERTIES
+            LINK_FLAGS "-undefined dynamic_lookup")
+    ELSE()
+        TARGET_LINK_LIBRARIES(_${name} ${PYTHON_LIBRARIES})
+    ENDIF()
 
     IF (NEKPY_LIBDEPENDS)
         TARGET_LINK_LIBRARIES(_${name} ${NEKPY_LIBDEPENDS})

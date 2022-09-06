class Nektar < Formula
  desc "High-performance spectral/hp element framework"
  homepage "https://www.nektar.info/"
  url "https://gitlab.nektar.info/nektar/nektar/-/archive/v5.2.0/nektar-v5.2.0.tar.bz2"
  sha256 "b58f7cff1d2579822c3f11a9b2bb0faa25c15709a8033755239a872ab5889719"

  bottle do
    root_url "https://github.com/mdave/homebrew-nektar/releases/download/nektar-5.2.0"
    sha256 cellar: :any, monterey: "08572273b055ceef537a161114ad5e60876eda9f4d24f8a0088602168a877aee"
    sha256 cellar: :any, big_sur:  "93ff190edd12c1b91e440a1daff98291d0f1826fe1212961921483ca6d49c1e5"
    sha256 cellar: :any, catalina: "723dadc6965eb8b97e9aeea40ecfafdd59e03f4c2cbb3ed7a88a5c56c5afa729"
  end

  depends_on "arpack"
  depends_on "boost"
  depends_on "boost-python3"
  depends_on "cmake"
  depends_on "fftw"
  depends_on "hdf5-mpi"
  depends_on "numpy"
  depends_on "open-mpi"
  depends_on "opencascade"
  depends_on "python@3.10"
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
    args << "-DPYTHON_EXECUTABLE=#{Formula["python@3.10"].opt_bin}/python3"

    # Compile with C++14 support for boost 1.75+
    args << "-DCMAKE_CXX_STANDARD=14"

    mkdir "build" do
      system "cmake", "..", *args
      system "make", "install"

      # Also need to install NekPy bindings
      python = Formula["python@3.10"]
      system python.bin/"python3", *Language::Python.setup_install_args(prefix)

      site_packages = Language::Python.site_packages(python)
      pth_contents = "import site; site.addsitedir('#{libexec/site_packages}')\n"
      (prefix/site_packages/"homebrew-nektar.pth").write pth_contents
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
      set(CMAKE_CXX_STANDARD 14)
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
    python = Formula["python@3.10"]
    system python.bin/"python3", "project.py"
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
diff --git a/cmake/ThirdPartyHDF5.cmake b/cmake/ThirdPartyHDF5.cmake
index c062dd833..a1983abf4 100644
--- a/cmake/ThirdPartyHDF5.cmake
+++ b/cmake/ThirdPartyHDF5.cmake
@@ -69,6 +69,10 @@ IF (NEKTAR_USE_HDF5)
         MESSAGE(STATUS "Found HDF5: ${HDF5_LIBRARIES}")
         SET(HDF5_CONFIG_INCLUDE_DIR ${HDF5_INCLUDE_DIRS})
         ADD_CUSTOM_TARGET(hdf5-1.8.16 ALL)
+
+        IF(HDF5_VERSION VERSION_GREATER_EQUAL 1.10.0)
+            ADD_DEFINITIONS(-DH5_USE_110_API)
+        ENDIF()
     ENDIF()
 
     MARK_AS_ADVANCED(HDF5_LIBRARIES)
diff --git a/library/NekMesh/CADSystem/OCE/OpenCascade.h b/library/NekMesh/CADSystem/OCE/OpenCascade.h
index 8ed19f9a2..47df00dd7 100644
--- a/library/NekMesh/CADSystem/OCE/OpenCascade.h
+++ b/library/NekMesh/CADSystem/OCE/OpenCascade.h
@@ -60,7 +60,6 @@
 #include <BRep_Tool.hxx>
 #include <GCPnts_AbscissaPoint.hxx>
 #include <GProp_GProps.hxx>
-#include <GeomAdaptor_HSurface.hxx>
 #include <GeomLProp_CLProps.hxx>
 #include <GeomLProp_SLProps.hxx>
 #include <ShapeAnalysis_Curve.hxx>
diff --git a/library/SolverUtils/Advection/Advection.h b/library/SolverUtils/Advection/Advection.h
index 277008d4e..3444b2c09 100644
--- a/library/SolverUtils/Advection/Advection.h
+++ b/library/SolverUtils/Advection/Advection.h
@@ -45,8 +45,6 @@
 #include <SolverUtils/SolverUtilsDeclspec.h>
 #include <iomanip>
 
-using namespace std;
-
 namespace Nektar
 {
 namespace SolverUtils
diff --git a/solvers/IncNavierStokesSolver/EquationSystems/StandardExtrapolate.cpp b/solvers/IncNavierStokesSolver/EquationSystems/StandardExtrapolate.cpp
index 6b75d7c71..ab1d2c8de 100644
--- a/solvers/IncNavierStokesSolver/EquationSystems/StandardExtrapolate.cpp
+++ b/solvers/IncNavierStokesSolver/EquationSystems/StandardExtrapolate.cpp
@@ -167,7 +167,7 @@ void StandardExtrapolate::v_AccelerationBDF(
         Array<OneD, NekDouble> accelerationTerm(nPts, 0.0);
         if (m_pressureCalls > 2)
         {
-            int acc_order = min(m_pressureCalls - 2, m_intSteps);
+            int acc_order = std::min(m_pressureCalls - 2, m_intSteps);
             Vmath::Smul(nPts, DuDt_Coeffs[acc_order - 1][0], array[0], 1,
                         accelerationTerm, 1);
 
diff --git a/library/FieldUtils/ProcessModules/ProcessInterpField.cpp b/library/FieldUtils/ProcessModules/ProcessInterpField.cpp
index b1cf6ec78..eab7a21d4 100644
--- a/library/FieldUtils/ProcessModules/ProcessInterpField.cpp
+++ b/library/FieldUtils/ProcessModules/ProcessInterpField.cpp
@@ -33,7 +33,6 @@
 ////////////////////////////////////////////////////////////////////////////////
 #include <iostream>
 #include <string>
-using namespace std;
 
 #include <boost/core/ignore_unused.hpp>
 #include <boost/geometry.hpp>
@@ -46,6 +45,7 @@ using namespace std;
 
 #include "ProcessInterpField.h"
 
+using namespace std;
 namespace bg  = boost::geometry;
 namespace bgi = boost::geometry::index;
 

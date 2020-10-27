//----------------------------------------------------------------------------------------------------------------------
// This declarative Jenkins pipeline encodes all the steps required for the nightly/continuous of a single platform.
// Other jobs may call this pipeline to execute the build, test and installation of a set platforms.
//
// Author: Pere Mato
//----------------------------------------------------------------------------------------------------------------------

pipeline {
  parameters {
    string(name: 'EXTERNALS', defaultValue: 'devgeantv/latest', description: 'LCG software stack in CVMFS')
    choice(name: 'MODE', choices: ['experimental', 'nightly'], description: 'CDash mode')
    string(name: 'ExtraCMakeOptions', defaultValue: '', description: 'CMake extra configuration options')
    string(name: 'VERSION', defaultValue: 'master', description: 'Branch to use for the build (master, experimental, ... )')
    choice(name: 'LABEL', choices: ['centos7', 'centos8', 'slc6', 'mac1015', 'cuda10', 'ubuntu18', 'ubuntu20'])
    choice(name: 'COMPILER', choices: ['gcc7', 'gcc8', 'gcc9', 'gcc10', 'clang8', 'clang10', 'native'])
    choice(name: 'BUILDTYPE', choices: ['Release', 'Debug'])
    choice(name: 'OPTION', choices: ['default', 'SPEC', 'AVX', 'GDML'])
    choice(name: 'BACKEND', choices: ['scalar', 'vc'])
    string(name: 'DOCKER_LABEL', defaultValue: 'docker-host', description: 'Label for the the nodes able to launch docker images')
  }

  environment {
    CMAKE_INSTALL_PREFIX = 'install'
    CMAKE_SOURCE_DIR     = 'vecgeom'
    CMAKE_BINARY_DIR     = 'build'
  }

  agent none

  stages {
    //------------------------------------------------------------------------------------------------------------------
    //---Build & Test stages--------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------
    stage('Prepa'){
      steps {
        init()
      }
    }
    stage('InDocker') {
      when {
        beforeAgent true
        expression { params.LABEL =~ 'centos|ubuntu' }
      }
      agent {
        docker {
          image "gitlab-registry.cern.ch/sft/docker/$LABEL"
          label "$DOCKER_LABEL"
          args  """-v /cvmfs:/cvmfs 
                   -v /ccache:/ccache 
                   -v /ec:/ec
                   -e SHELL
                   --net=host
                   --hostname ${LABEL}-docker
                """
        }
      }
      stages {
        stage('Build&Test') {
          steps {
            buildAndTest()
          }
          post {
            success {
              deleteDir()
            }
          }
        }
      }
    }
    stage('InGPU') {
      when {
        beforeAgent true
        expression { params.LABEL =~ 'cuda' }
      }
      agent {
        label 'cuda10'
      }
      stages {
        stage('Build&Test') {
          steps {
            buildAndTest()
          }
          post {
            success {
              deleteDir()
            }
          }
        }
      }
    }
  }
}

def init() {
  currentBuild.displayName = "#${BUILD_NUMBER}" + ' ' + params.OPTION + '-' + params.BACKEND + '-' + params.LABEL + '-' + params.COMPILER + '-' + params.BUILDTYPE
}

def buildAndTest() {
  sh label: 'build_and_test', script: """
    source /cvmfs/sft.cern.ch/lcg/views/${EXTERNALS}/x86_64-centos7-${COMPILER}-opt/setup.sh
    env | sort | sed 's/:/:?     /g' | tr '?' '\n'
    ctest -V -S vecgeom/jenkins/vecgeom-ctest.cmake,$MODE
  """
}
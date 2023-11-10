pipeline {
    agent { label 'jenkinsfile' }
    triggers {
        pollSCM('H/10 * * * *')
    }

    options {
        disableConcurrentBuilds()
        buildDiscarder(logRotator(numToKeepStr: '50'))
        timestamps()
    }

    stages {
        stage ('Chameleon only') {
            steps {
                sh '''#!/bin/bash -le
                    ####################################################
                    # Configure and build
                    ####################################################
                    module purge
                    module load gcc/10.2.0
                    module load cmake/3.21.2
                    ####################################################
                    # BLAS/LAPACK
                    ####################################################
                    module load mkl/2020.0.166
                    ####################################################

                    set -x
                    ./config.sh -t -e
                    ./clean_build.sh
                '''
            }
        }
        stage ('test Chameleon') {
            steps {

                sh '''#!/bin/bash -le
                    ####################################################
                    # Run tester
                    ####################################################
                    echo "========================================"
                    module purge
                    module load gcc/10.2.0
                    module load cmake/3.21.2
                    ####################################################
                    # BLAS/LAPACK
                    ####################################################
                    module load mkl/2020.0.166
                    cd bin/
                    ctest --no-compress-output --verbose
                    '''
            }
        }
        stage ('HiCMA') {
            steps {
                sh '''#!/bin/bash -le
                    ####################################################
                    # Configure and build
                    ####################################################
                    module purge
                    module load gcc/10.2.0
                    module load cmake/3.21.2
                    ####################################################
                    # BLAS/LAPACK
                    ####################################################
                    module load mkl/2020.0.166
                    ####################################################

                    set -x
                    ./config.sh -t -e -H
                    ./clean_build.sh
                '''
            }
        }
        stage ('test HiCMA') {
            steps {

                sh '''#!/bin/bash -le
                    ####################################################
                    # Run tester
                    ####################################################
                    echo "========================================"
                    module purge
                    module load gcc/10.2.0
                    module load cmake/3.21.2
                    ####################################################
                    # BLAS/LAPACK
                    ####################################################
                    module load mkl/2020.0.166
                    module load gsl/2.6-gcc-10.2.0
                    cd bin/
                    ctest --no-compress-output --verbose
                    '''
            }
        }
	    stage('documentation') {
             agent { label 'jenkinsfile'}
             steps {
                 sh '''#!/bin/bash -le
                    module purge
                    module load gcc/10.2.0
                    module load cmake/3.21.2
                    ####################################################
                    # BLAS/LAPACK
                    ####################################################
                    module load mkl/2020.0.166
                    module load gsl/2.6-gcc-10.2.0
                    ./config.sh -t -e
                    ./clean_build.sh
                    cd bin
                    make docs
                    '''
                 publishHTML( target: [allowMissing: false, alwaysLinkToLastBuild: false, keepAll: false, reportDir: 'bin/docs/bin/html', reportFiles: 'index.html', reportName: 'Doxygen Documentation', reportTitles: ''] )
             }
        }
    }

    // Post build actions
    post {
        //always {
        //}
        //success {
        //}
        //unstable {
        //}
        //failure {
        //}
        unstable {
                emailext body: "${env.JOB_NAME} - Please go to ${env.BUILD_URL}", subject: "Jenkins Pipeline build is UNSTABLE", recipientProviders: [[$class: 'CulpritsRecipientProvider'], [$class: 'RequesterRecipientProvider']]
        }
        failure {
                emailext body: "${env.JOB_NAME} - Please go to ${env.BUILD_URL}", subject: "Jenkins Pipeline build FAILED", recipientProviders: [[$class: 'CulpritsRecipientProvider'], [$class: 'RequesterRecipientProvider']]
        }
    }
}
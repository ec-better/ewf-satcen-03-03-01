<<<<<<< HEAD
peline {
=======
node('ci-community') {
  
  stage 'Checkout'
  checkout scm
  
  stage 'Setup environment'
  env.PATH = "${tool 'apache-maven-3.0.5'}/bin:/opt/anaconda/bin:${env.PATH}"
  
  stage 'Package and Deploy'
  sh 'mvn deploy -Drelease=true'
>>>>>>> master

  options {
    buildDiscarder(logRotator(numToKeepStr: '5'))
  }

  environment {
        PATH="/opt/anaconda/bin:$PATH"
  }

  agent {
    node {
      label 'ci-community-docker'
    }
  }

  stages {

    stage('Package & Dockerize') {
      steps {
        withMaven( maven: 'apache-maven-3.0.5' ) {
            sh 'mvn -B deploy'
        }
      }
    }
  }
}


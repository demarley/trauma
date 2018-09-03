#ifndef DEEPLEARNING_H
#define DEEPLEARNING_H

#include <string>
#include <map>
#include <vector>

#include "lwtnn/lwtnn/interface/LightweightNeuralNetwork.hh"
#include "lwtnn/lwtnn/interface/parse_json.hh"

#include "trauma/miniAna/interface/tools.h"
#include "trauma/miniAna/interface/configuration.h"
#include "trauma/miniAna/interface/physicsObjects.h"


class deepLearning {
  public:
    deepLearning( configuration& cmaConfig );

    ~deepLearning();

    void training(const Track& tk, const Muon& mu);
    void inference(const Track& tk, const Muon& mu);
    void loadFeatures();

    std::map<std::string,double> predictions() const {return m_discriminant;}
    double prediction();
    double prediction(const std::string& key);
    std::map<std::string,double> features() const {return m_features;}

  protected:

    configuration *m_config;

    lwt::LightweightNeuralNetwork* m_lwnn;       // LWTNN tool
    std::map<std::string, double> m_features;    // values for inputs to the DNN
    std::string m_dnnKey;                        // default key for accessing map of values
    float m_DNN;                                 // DNN prediction for one key

    Track m_track;
    Muon m_muon;

    std::map<std::string,double> m_discriminant; // map of DNN predictions
};

#endif

//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <string>
#include <memory>
#include <vector>
#include <limits>
#include "../thirdParty/cget/include/savvy/reader.hpp"
#include "../thirdParty/cget/include/savvy/armadillo_vector.hpp"

using namespace std;


class Genotypes {
protected:
   string file;
   unique_ptr<savvy::vcf::indexed_reader<1>> f;
   bool is_reading;
   bool has_cached;
   savvy::site_info anno;
   savvy::armadillo::dense_vector<double> dosages;
   bool flip;

public:
   Genotypes(const Rcpp::CharacterVector& v ): file(v[0]), is_reading(false), has_cached(false), flip(false) {}

   void open(const vector<string>& samples) {
      f = unique_ptr<savvy::vcf::indexed_reader<1>>(new savvy::vcf::indexed_reader<1>(file, {""}, savvy::fmt::ds));
      if (samples.size() > 0) {
         f->subset_samples({ samples.begin(), samples.end() });
      }
      is_reading = false;
      has_cached = false;
   }

   Rcpp::CharacterVector get_samples() const {
      return Rcpp::wrap(savvy::vcf::reader<1>(file, savvy::fmt::gt).samples());
   }

   Rcpp::CharacterVector get_chromosomes() const {
      return Rcpp::wrap(savvy::vcf::indexed_reader<1>(file, {""}, savvy::fmt::gt).chromosomes());
   }

   arma::vec read_variant(const string& chromosome, unsigned int position, const string& risk_allele, const string& protective_allele) {
      if (is_reading) {
         if ((chromosome.compare(anno.chromosome()) != 0) || (position > anno.position() + 10000)) {
            f->reset_region({ chromosome, position, numeric_limits<int>::max() });
            has_cached = false;
         }
      } else {
         f->reset_region({ chromosome, position, numeric_limits<int>::max() });
         is_reading = true;
         has_cached = false;
      }
      if (!has_cached) {
         f->read(anno, dosages);
      }
      bool position_found = false;
      while (f->good()) {
         if (anno.position() < position) {
            f->read(anno, dosages);  
            continue;
         }
         if (anno.position() > position) {
            has_cached = true;
            break;
         }
         position_found = true;
         flip = false;
         if ((risk_allele.compare(anno.ref()) == 0) && (protective_allele.compare(anno.alt()) == 0)) {
            flip = true;
         } else if ((risk_allele.compare(anno.alt()) != 0) or (protective_allele.compare(anno.ref()) != 0)) {
            f->read(anno, dosages);
            continue;
         }
         if (flip) {
            dosages -= 2.0;
            dosages *= -1;
         }
         has_cached = false;
         return dosages;
      }
      if (position_found) {
         Rcpp::warning("Allele combination %s/%s not found at position %s:%d", protective_allele, risk_allele, anno.chromosome(), anno.position());
      } /* else {
         Rcpp::warning("No variant at position %s:%d", anno.chromosome(), anno.position());
      } */
      return arma::vec(); // return empty vector if variant was not found
   }
};


RCPP_MODULE(GenotypesEx) {
   Rcpp::class_<Genotypes>("Genotypes")
      .constructor<Rcpp::CharacterVector>()
      .method("open", &Genotypes::open)
      .method("get_samples", &Genotypes::get_samples)
      .method("get_chromosomes", &Genotypes::get_chromosomes)
      .method("read_variant", &Genotypes::read_variant)
      ;
}

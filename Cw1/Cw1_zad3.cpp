/*! Laboratorium: Zaawansowane przetwarzanie sygnałów telekomunikacji cyfrowej
 *  Ćw. 1. zad. 3.
 *   - wczytanie z pliku odpowiedzi impulsowej filtru
 *    (skrypt generujący plik danych: ./matlab/design_zad_2.m
 *   - klasyczny interpolator L = 2
 *
 * \author Marek Blok
 * \date 2018.03.01
 */
#include <DSP_lib.h>

int main(int argn, char *args[])
{
  /*************************************************************/
  // Log file setup
  DSP::log.SetLogFileName("log_file.txt");
  DSP::log.SetLogState(DSP::e::LogState::file | DSP::e::LogState::console);

  DSP::log << DSP::lib_version_string() << std::endl << std::endl;
  /*************************************************************/

  DSP::LoadCoef coef_info;
  int N_LPF;
  DSP::Float_vector h_LPF;

  // wczytanie wspó³czynników z pliku matlab/cw1_zad2.coef
  coef_info.Open("cw1_zad2.coef", "matlab");
  N_LPF = coef_info.GetSize(0); // odczytanie liczby wspó³czynników pierwszeego wektora (indeksowanie od zera)
  if (N_LPF < 1)
  {
    DSP::log << DSP::e::LogMode::Error << "No filter coefficients: aboarding" << std::endl;
    return -1;
  }
  else
  {
    // wczytanie N_LPF wspó³czynników do tablicy h_LPF
    coef_info.Load(h_LPF);
  }
  /*************************************************************/

  DSP::Clock_ptr MasterClock;
  // Pozyskanie wskazania do g³ównego zegara algorytmu
  MasterClock=DSP::Clock::CreateMasterClock();

  long int Fp1, Fp2;

  // Ÿród³o próbek: plik ./test.wav (bie¿¹cy katalog)
  DSP::u::WaveInput AudioIn(MasterClock, "test.wav", ".");
  // odczytanie szybkoœci próbkowania pliku powi¹zanego z bloczkiem AudioIn
  Fp1 = AudioIn.GetSamplingRate();

  Fp2 = 2*Fp1;
  // dwukrotny zeroinserter z wejœciem pracuj¹cym w rytm zegara MasterClock
  DSP::u::Zeroinserter Zeroinserter(MasterClock, 2U);
  // przyk³ad dodania tekstu do domyslnej nazwy bloczka (dla pliku *.dot)
  Zeroinserter.SetName("x2");

  // filtr typu FIR o rzeczywistych wspó³czynnikach odpowiedzi impulsowej (tutaj filtr interpolacyjny)
  // d³ugoœæ odpowiedzi impulsowej N_LPF, próbki odpowiedzi impulsowej w tablicy h_LPF
  DSP::u::FIR InterpFIR(h_LPF);

  // wyjœcie na kartê dŸwiêkow¹ - szybkoœæ próbkowania Fp2 (mono/16bit)
  DSP::u::AudioOutput SoundOut(Fp2, 1, 16);
  // zapis do pliku *.wav (mono/16bit, szybkosæ próbkowania Fp2)
  DSP::u::FileOutput FileOut("cw1_zad3.wav", DSP::e::SampleType::ST_short, 1, DSP::e::FileType::FT_wav, Fp2);

  /*************************************************************/
  // Definicje po³¹czeñ pomiêdzy bloczkami
  // wyjœcie "out" bloczka AudioIn ³¹czymy w wejœciem "in" bloczka Zeroinserter
  AudioIn.Output("out") >> Zeroinserter.Input("in");
  // wyjœcie "out" bloczka Zeroinserter ³¹czymy w wejœciem "in" bloczka InterpFIR
  Zeroinserter.Output("out") >> InterpFIR.Input("in");
  // wyjœcie "out" bloczka InterpFIR ³¹czymy w wejœciem "in" bloczka SoundOut
  InterpFIR.Output("out") >> SoundOut.Input("in");
  // wyjœcie "out" bloczka InterpFIR ³¹czymy w wejœciem "in" bloczka FileOut
  // Uwaga: jedno wyjœcie mo¿na pod³¹czyæ do kilku wejœæ ale nie mo¿na
  //   pod³¹czyæ kilku wyjœæ do jednego wejœcia
  InterpFIR.Output("out") >> FileOut.Input("in");


  /////////////////////////////////
  // sprawdzenie, czy do wszystkich wejsæ, wszystkich bloczków pod³¹czone s¹ sygna³y wejsciowe
  DSP::Component::CheckInputsOfAllComponents();

  // *********************************** //
  // Zapis schematu zaimplementowanego algorytmu do pliku *.dot
  // (przetworzenie do pliku *.gif skryptem rundot.bat)
  DSP::Clock::SchemeToDOTfile(MasterClock, "Cw1_zad3.dot");

  // *********************************** //
  // przetwarzany w blokach po SamplesInSegment próbek wejsciowych
  int SamplesInSegment = 512;

  // licznik przetworzonych próbek wejsciowych
  __int64 NoOfSamplesProcessed = 0;

  // Maksymalna liczba przetwarzanych próbek wejsciowych (10 sekund)
  #define MAX_SAMPLES_TO_PROCESS 10*Fp1

  while(NoOfSamplesProcessed < MAX_SAMPLES_TO_PROCESS)
  {

    // ********************************************************** //
    // uruchom zaimplementowany algorytm na SamplesInSegment cykli zegara MasterClock
    DSP::Clock::Execute(MasterClock, SamplesInSegment);
    // ********************************************************** //

    if (AudioIn.GetBytesRead() > 0)
    { // resetuj licznik przetworzonych próbek je¿eli uda³o siê odczytaæ próbki z pliku
      // - wymusza przetworzenie ca³ego pliku wejsciowego
      NoOfSamplesProcessed = 0;
    }
    else // Play 200ms more
    { // je¿eli nie odczytano próbek z pliku wejsciowego to przetwórz jeszcze tylko 1/5 sekundy
      if (NoOfSamplesProcessed < MAX_SAMPLES_TO_PROCESS - Fp1/5)
          NoOfSamplesProcessed = MAX_SAMPLES_TO_PROCESS - Fp1/5;
    }

    NoOfSamplesProcessed += SamplesInSegment;
    // ********************************************************** //
  }

  /*************************************************************/
  // Zwolnienie zarezerwowanych zasobów
  DSP::Clock::FreeClocks();
  /*************************************************************/

  return 0;
}

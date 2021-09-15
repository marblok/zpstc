/*! Laboratorium: Zaawansowane przetwarzanie sygna��w telekomunikacji cyfrowej
 *  �w. 6. zad. 2.
 *
 * \author Marek Blok
 * \date 2018.03.16
 */

#include <sstream>
using namespace std;

#include <DSP_lib.h>
#include <memory>

DSP::Float_ptr read_buffer = NULL;
int buffer_size;
int No_of_samples=0;
// //int K = 16; int M = 4; // Oversampling ratio
//int K = 32; int M = 8; // Oversampling ratio : K/M
int K = 32; int M = 1; // Oversampling ratio : K/M


int main(int argn, char *args[])
{
    /*! Koncepcja
    * - cztery kana�y - r�ne modulacje - ten sam filtr kszta�tuj�cy
    *   DPSK, pi/4-QPSK, (MSK), OQPSK, QAM-16, OOK
    * - FFT - 16 ==> od razu lokowanie na PCz (nadpor�bkowanie 16)
    * - dost�pne kana�y 1..7 (nie mo�na u�y� dla FMT dw�ch s�siednich ze wzgl�du na ICI)
    *   u�ycie: 2, 5, 7
    * reszta sta�e zero !!! SetConstInput (chyba powinno dzia�a�)
    */
  /*************************************************************/
  // Log file setup
  DSP::log.SetLogFileName("log_file.txt");
  DSP::log.SetLogState(DSP::e::LogState::file | DSP::e::LogState::console);

  DSP::log << DSP::lib_version_string() << std::endl << std::endl;
  /*************************************************************/

  long int Fp2, F_symb;
  DSP::LoadCoef coef_info;
  int N_rc;
  //unsigned int L;
  DSP::Float_vector h_rc;
  string fft_output_name;
  string tekst, tekst2;

  map<string,shared_ptr<DSP::Component> > blocks;

  coef_info.Open("ex6_task1_h_rc.coef", "matlab");
  N_rc = coef_info.GetSize(0);
  if (N_rc < 1)
  {
    DSP::log << DSP::e::LogMode::Error << "No ex6_task1_h_rc.coef: aboarding" << std::endl;
    return -1;
  }
  else
  {
    coef_info.Load(h_rc);
    for (auto n=0U; n<h_rc.size(); n++){
        h_rc[n] *= K;
    }
    F_symb = coef_info.Fp;
  }
  /*************************************************************/

  DSP::Clock_ptr SymbolClock, InputClock;
  InputClock=DSP::Clock::CreateMasterClock();


  //DSP::u::WaveInput AudioIn(MasterClock, "test.wav", ".");
  //F_symb = AudioIn.GetSamplingRate();


  F_symb = 1500; // dla L = 32 => Fp2 = 48000
  Fp2 = K*F_symb;

  DSP::log << "Fsymb = " << int(F_symb) << ", Fp2 = " << int(Fp2) << ", K = " << K << ", M = " << M << std::endl;

  SymbolClock=DSP::Clock::GetClock(InputClock, 1, K/M);


  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
  // Multichannel input signal (for complex valued signal (2U))
  string input_filename;

  char mode = 'A';
  if (argn == 2) {
    mode = args[1][0];
  }
  unsigned int input_delay = 1U;
  // 1 symbol delay at on the transmitting side - because of FFT
  // 2 bit delay after decoder
  // previous: SetSkip for 10 dibits
  // new: 10 dibits + 1 dibit (additional symbol from Serial2Parallel converter located before DSP::u::SymbolMapper)
  unsigned int bin_data_skip = (2*(1+ 2+8));

  switch (mode) {
    case 'A':
        input_filename = "ex6_task1.flt";
        break;
    case 'B':
        input_filename = "ex6_task1b.flt";
        // skip additional symbol (compensate for the DSP::u::Parallel2Serial in the modulator)
        input_delay += K;
        bin_data_skip += 2*1; // skip additional dibit
        break;
    default:
        DSP::log << DSP::e::LogMode::Error << "Unsupported mode" << std::endl;
        return -1;
        break;
  }
  blocks["InputSignal"] = make_shared<DSP::u::FileInput>(InputClock, input_filename, 2U, DSP::e::SampleType::ST_float, DSP::e::FileType::FT_flt);
  blocks["InputSignal"]->SetName(input_filename);
  blocks["DumpImag"] = make_shared<DSP::u::Vacuum>();

  // "Frame"/ FFT symbol synchronization 
  blocks["SymbolTimigDelay"] = shared_ptr<DSP::Block>(new DSP::u::Delay (input_delay));

  blocks["InputSignal"]->Output("out.re") >> blocks["SymbolTimigDelay"]->Input("in");
  blocks["InputSignal"]->Output("out.im") >> blocks["DumpImag"]->Input("in");

  blocks["OutputBuffer"] = make_shared<DSP::u::Serial2Parallel>(InputClock, K,1);
  blocks["SymbolTimigDelay"]->Output("out") >> blocks["OutputBuffer"]->Input("in");

  blocks["fft"] = make_shared<DSP::u::FFT>(K, false);
  // poliphase filters
  //DSP::u::FIR *H_g[K];
  vector<shared_ptr<DSP::u::FIR> > H_g(K, nullptr);

  for (int ind = 0; ind < K; ind++)
  {
    H_g[ind] = make_shared<DSP::u::FIR>(h_rc, (K-1)-ind, K, M);

    tekst = "out" + to_string(ind+1);
    tekst2 = "in" + to_string(ind+1);
    blocks["OutputBuffer"]->Output(tekst) >> H_g[ind]->Input("in");
    H_g[ind]->Output("out") >> blocks["fft"]->Input(tekst2);
  }

  // number of active subchannels
  const int NoOfActiveChannels = 3;

  vector<shared_ptr<DSP::u::Vacuum> >         Discard(K-NoOfActiveChannels, nullptr);
  vector<shared_ptr<DSP::u::RawDecimator> >   RawDec(NoOfActiveChannels, nullptr);
  vector<shared_ptr<DSP::u::SymbolDemapper> >  PSKdecoder(NoOfActiveChannels, nullptr);
  vector<shared_ptr<DSP::u::Parallel2Serial> > PSconv(NoOfActiveChannels, nullptr);
  vector<shared_ptr<DSP::u::FileOutput> >     SymbData(NoOfActiveChannels, nullptr);
  vector<shared_ptr<DSP::u::FileOutput> >     BinData(NoOfActiveChannels, nullptr);

  // channel1 = 1; channel2 = 4; channel3 = 6;
  // subchannels numbering convension: 0, 1, ..., K-1
  int channel_no[NoOfActiveChannels] = {1, 4, 6};

  int current_discard_block_no = 0;
  int active_channel_index = -1;
  for (int ind = 0; ind < K; ind++)
  {
    active_channel_index = -1;
    for (int ind2 = 0; ind2 < NoOfActiveChannels; ind2++)
    {
      if (ind == channel_no[ind2])
      { // memorise indexes of consecative active channels (indexing from zero)
        active_channel_index = ind2;
      }
    }

    fft_output_name = "out" + to_string(ind+1);
    if (active_channel_index == -1)
    { // subchannel is not on the list of active channels
      Discard[current_discard_block_no] = make_shared<DSP::u::Vacuum>(true);
      //Discard[current_discard_block_no] = shared_ptr<DSP::u::Vacuum>(new DSP::u::Vacuum(true));
      blocks["fft"]->Output(fft_output_name) >> Discard[current_discard_block_no++]->Input("in");
    }
    else
    { // subchannel is on the list of active channels
      stringstream filename;

//      PSKdecoder[active_channel_index] = make_shared<DSP::u::SymbolDemapper>(SymbolClock, DSP::e::ModulationType::PSK, 2, 0.0); //M_PIx1f/2); // QPSK_A
      PSKdecoder[active_channel_index] = make_shared<DSP::u::SymbolDemapper>(DSP::e::ModulationType::PSK, 2, 0.0); // QPSK_A
      PSconv[active_channel_index] = make_shared<DSP::u::Parallel2Serial>(SymbolClock, 2);
      filename << "ex6_task2_symb_" << mode << "_" << char('a' + active_channel_index) << ".flt";
      SymbData[active_channel_index] = make_shared<DSP::u::FileOutput>(filename.str(), DSP::e::SampleType::ST_float, 2U, DSP::e::FileType::FT_flt, F_symb);
      SymbData[active_channel_index]->SetName(filename.str());

      filename.str(""); filename.clear();
      filename << "ex6_task2_bin_" << mode << "_" << char('a' + active_channel_index) << ".txt";
      BinData[active_channel_index] = make_shared<DSP::u::FileOutput>(filename.str(), DSP::e::SampleType::ST_bit, 1U, DSP::e::FileType::FT_raw);
      //BinData[active_channel_index] = make_shared<DSP::u::FileOutput>(filename.str(), DSP::e::SampleType::ST_float, 1U, DSP::e::FileType::FT_raw);
      BinData[active_channel_index]->SetName(filename.str());

      // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
      // dibit offset - byte start synchronization
      BinData[active_channel_index]->SetSkip(bin_data_skip);
      // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

      // FFT with poliphase filters (only for used channels)
      if (M > 1)
      {
        RawDec[active_channel_index] = make_shared<DSP::u::RawDecimator>(SymbolClock, M, 2U);
        blocks["fft"]->Output(fft_output_name) >>  RawDec[active_channel_index]->Input("in");
        RawDec[active_channel_index]->Output("out") >> SymbData[active_channel_index]->Input("in");
        RawDec[active_channel_index]->Output("out") >> PSKdecoder[active_channel_index]->Input("in");
      }
      else
      {
        RawDec[active_channel_index] = NULL;
        blocks["fft"]->Output(fft_output_name) >> SymbData[active_channel_index]->Input("in");
        blocks["fft"]->Output(fft_output_name) >> PSKdecoder[active_channel_index]->Input("in");
      }
      PSKdecoder[active_channel_index]->Output("out") >> PSconv[active_channel_index]->Input("in");
      PSconv[active_channel_index]->Output("out") >> BinData[active_channel_index]->Input("in");
    }
  }



  // *********************************** //
  DSP::Clock::SchemeToDOTfile(SymbolClock, "ex6_task2.dot");
  /////////////////////////////////
  // check if there are signals
  // connected to all inputs
  DSP::Component::CheckInputsOfAllComponents();

  // *********************************** //
  int SamplesInSegment = 512;
  __int64 NoOfSamplesProcessed = 0;
  // 10 seconds
  #define MAX_SAMPLES_TO_PROCESS 50000 //1*Fp2
  while(NoOfSamplesProcessed < MAX_SAMPLES_TO_PROCESS)
  {

    // ********************************************************** //
    DSP::Clock::Execute(InputClock, SamplesInSegment);
    // ********************************************************** //

    /*
    if (BinData1.GetBytesRead() > 0)
    {
        NoOfSamplesProcessed = 0; // Play the whole file
    }
    else // Play 200ms more
    {
      if (NoOfSamplesProcessed < MAX_SAMPLES_TO_PROCESS - Fp2/5)
          NoOfSamplesProcessed = MAX_SAMPLES_TO_PROCESS - Fp2/5;
    }
    */

    NoOfSamplesProcessed += SamplesInSegment;
    // ********************************************************** //

    //sprintf(tekst, "NoOfSamplesProcessed = %i", int(NoOfSamplesProcessed));
    //DSPf_InfoMessage(tekst);
 }

  /*************************************************************/
  H_g.clear();
  blocks.clear();
  Discard.clear();
  PSKdecoder.clear();
  RawDec.clear();
  SymbData.clear();
  BinData.clear();

  DSP::Clock::ListComponents();
  /*************************************************************/
  DSP::Clock::FreeClocks();

  return 0;
}

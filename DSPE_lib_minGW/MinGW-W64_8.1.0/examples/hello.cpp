/*! Simple Digital Signal Processing Engine usage example.
 * \author Marek Blok
 * \date 2006.09.11
 * \date updated 2021.01.18
 */
#include <DSP_lib.h>

#ifndef INCLUDE_DSPE_EXAMPLES
int main(void)
#else
#include "DSPE_examples.h"
int test_hello(void)
#endif // INCLUDE_DSPE_EXAMPLES
{
  DSP::Clock_ptr MasterClock;
  string tekst;
  int temp;
  long int Fp;

  DSP::log.SetLogState(DSP::e::LogState::console | DSP::e::LogState::file);
  DSP::log.SetLogFileName("log_file.log");

  DSP::log << DSP::lib_version_string() << endl << endl;
  DSP::log << "Hello" << DSP::e::LogMode::second << "World !!!" << endl;
 
  MasterClock=DSP::Clock::CreateMasterClock();

  DSP::u::WaveInput AudioIn(MasterClock, "DSPElib.wav", ".");
  Fp = AudioIn.GetSamplingRate();

  DSP::u::AudioOutput AudioOut(Fp);

  AudioIn.Output("out") >> AudioOut.Input("in");

  DSP::Component::CheckInputsOfAllComponents();
  DSP::Clock::SchemeToDOTfile(MasterClock, "hello.dot");

  temp=1;
  do
  {
    DSP::Clock::Execute(MasterClock, Fp/8);

    DSP::log << "MAIN" << DSP::e::LogMode::second << temp << endl;
    temp++;
  }
  while (AudioIn.GetBytesRead() != 0);

  DSP::Clock::FreeClocks();
  DSP::log << DSP::e::LogMode::Error << "MAIN" << DSP::e::LogMode::second << "end" << endl;

  return 0;
}

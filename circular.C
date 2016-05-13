#include "circular.h"

// Function definitions
circ_buffer::circ_buffer(int size)
{
  idx = 0;
  buffer.resize(size);
}

void circ_buffer::push(std::string str) {
  buffer[idx] = str;
  idx = (idx + 1) % buffer.size();
}

// doesn't remove entry, just moves pointer to the next one
std::string circ_buffer::pop() {
  std::string r = buffer[idx];
  idx = (idx + 1) % buffer.size();
  return r;
}

void circ_buffer::flushToDisk(std::ofstream& logfile) {
  int i = idx;

  do {
    logfile << buffer[i];
    i = (i + 1) % buffer.size();
  } while (i != idx);

  logfile.flush();
}

void circ_buffer::pup(PUP::er& p) {
  p|buffer;
  p|idx;
}

void circ_buffer::resize(int size) {
  // Don't care about maintaining idx pointer correctness 
  // between resizes.
  // std::cout << "resizing from " << buffer.size() << " to " << size << std::endl;
  buffer.resize(size);
  idx = 0;
}


// void RingStream::flush() { #<{(| ignore regular flushes |)}># }
void RingStream::really_flush(std::ofstream& logfile) {
  buffer.rbuf.flushToDisk(logfile);
}

void RingStream::pup(PUP::er& p) {
  p|buffer;
}

void RingStream::open(const char *filename, ios_base::openmode mode) {}

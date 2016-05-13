#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <pup_stl.h>

#include <iostream>
#include <sstream>

class circ_buffer {
public:
  // Function definitions
  circ_buffer(int size = 100);

  void push(std::string str);

  // doesn't remove entry, just moves pointer to the next one
  std::string pop();

  void flushToDisk(std::ofstream& logfile);

  void pup(PUP::er& p);

  void resize(int size);

private:
  std::vector<std::string> buffer;
  int idx;
};

// Modified from:
//   http://stackoverflow.com/questions/2212776/overload-handling-of-stdendl
class RingStream: public std::ostream
{
public:
  RingStream(std::ostream& str)
    :std::ostream(&buffer)
    ,buffer(str)
  { }

  RingStream(int size = 100)
    :std::ostream(&buffer),
     buffer(foo)
  {
    // std::cout << "constructor called with size = " << size << std::endl;
    buffer.rbuf.resize(size);
  }

  // void flush() { #<{(| ignore regular flushes |)}># }

  // void RingStream::flush() { #<{(| ignore regular flushes |)}># }
  void really_flush(std::ofstream& logfile);

  void pup(PUP::er& p);

  void open(const char *filename, ios_base::openmode mode = ios_base::out);

private:
  // Write a stream buffer that prefixes each line with Plop
  class RingStreamBuf: public std::stringbuf
  {
    public:
      RingStreamBuf(std::ostream& str)
          :output(str)
      {}

      // When we sync the stream with the output. 
      // 1) Output Plop then the buffer
      // 2) Reset the buffer
      // 3) flush the actual output stream we are using.
      virtual int sync ( )
      {
        // std::cout << "sync called with str " << str() << std::endl;
        rbuf.push(str());
        str("");
        output.flush();
        return 0;
      }

      void pup(PUP::er& p) {
        p|rbuf;
        // p|output;
      }

      circ_buffer rbuf;
    private:
      std::ostream& output;
  };

  // My Stream just uses a version of my special buffer
  RingStreamBuf buffer;

  std::ostringstream foo;
};

///////////////////////////////////////////////////////////////////////
// Configuration file reader
// Author: Jeff Comer <jcomer2@illinois.edu>
#ifndef READER_H
#define READER_H

#include "useful.H"

class Reader {
public:
  Reader(const char* fileName) {
    FILE* inp = fopen(fileName, "r");
    char line[256];

    const int numLines = countParameterLines(fileName);
    param = new String[numLines];
    value = new String[numLines];

    int count = 0;
    while (fgets(line, 256, inp) != NULL) {
      // Ignore comments.
      int len = strlen(line);
      if (line[0] == '#') continue;
      if (len < 2) continue;
      
      String s(line);
      int numTokens = s.tokenCount();
      if (numTokens != 2) {
	printf("Warning: Invalid config file line: %s\n", line);
	continue;
      }
      
      String* tokenList = new String[numTokens];
      s.tokenize(tokenList);
      if (tokenList == NULL) {
	printf("Warning: Invalid config file line: %s\n", line);
	continue;
      }
      param[count] = tokenList[0];
      value[count] = tokenList[1];
      //printf("%s %s\n", tokenList[0].val(), tokenList[1].val());
      count++;

      delete[] tokenList;
    }
    num = count;

    fclose(inp); 
  }

  ~Reader() {
    delete[] param;
    delete[] value;
  }

  static int countParameterLines(const char* fileName) {
    FILE* inp = fopen(fileName, "r");
    char line[256];
    int count = 0;

    while (fgets(line, 256, inp) != NULL) {
      // Ignore comments.
      int len = strlen(line);
      if (line[0] == '#') continue;
      if (len < 2) continue;
      
      count++;
    }
    fclose(inp);

    return count;
  }

  int length() const { return num; }

  String getParameter(int i) const {
    i %= num;
    if (i < 0) i += num;

    return param[i];
  }

  String getValue(int i) const {
    i %= num;
    if (i < 0) i += num;

    return value[i];
  }
  
  String toString() const {
    String ret;

    for (int i = 0; i < num; i++) {
      ret.add(param[i]);
      ret.add(' ');
      ret.add(value[i]);
      ret.add('\n');
    }

    return ret;
  }

  int countParameter(const String& p) const {
    int count = 0;
    for (int i = 0; i < num; i++) {
      if (param[i] == p) count++;
    }
    return count;
  }

private:
  int num;
  String* param;
  String* value;
};

#endif

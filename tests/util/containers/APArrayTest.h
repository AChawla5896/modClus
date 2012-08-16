#ifndef AP_ARRAY_TEST_H
#define AP_ARRAY_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/containers/DArray.h>
#include <util/containers/APArray.h>
#include <util/containers/PArrayIterator.h>

using namespace Util;

class APArrayTest : public UnitTest 
{

private:

   const static int capacity = 10;

   typedef int Data;

   DArray<Data> array; 
   APArray<Data> parray;
   PArrayIterator<Data> iterator;
   
public:

   void setUp();
   void tearDown();
  
   void testAppend();
   void testAppendEmpty();
   void testModify();
   void testIterator();

};

void APArrayTest::setUp()
{
   array.allocate(capacity);
   for (int i=0; i < capacity; i++) {
      array[i] = (i+1)*10 + 1;
   }
}

void APArrayTest::tearDown()
{}
  
void APArrayTest::testAppend()
{
   printMethod(TEST_FUNC);

   parray.reserve(2);
   parray.append(array[8]); // 0
   parray.append(array[4]); // 1
   TEST_ASSERT(parray.capacity() == 2);
   parray.append(array[3]); // 2
   TEST_ASSERT(parray.capacity() == 4);
   parray.append(array[5]); // 3
   TEST_ASSERT(parray.capacity() == 4);
   parray.append(array[6]); // 4
   TEST_ASSERT(parray.capacity() == 8);
   parray.append(array[9]); // 5

   TEST_ASSERT(parray.size() == 6);
   TEST_ASSERT(parray[0] == array[8]);
   TEST_ASSERT(parray[1] == array[4]);
   TEST_ASSERT(parray[2] == array[3]);
   TEST_ASSERT(parray[3] == array[5]);
   TEST_ASSERT(parray[4] == array[6]);
   TEST_ASSERT(parray[5] == array[9]);

   parray.clear();
   TEST_ASSERT(parray.size() == 0);
   TEST_ASSERT(parray.capacity() == 8);

} 

void APArrayTest::testAppendEmpty()
{
   printMethod(TEST_FUNC);

   TEST_ASSERT(parray.size() == 0);
   TEST_ASSERT(parray.capacity() == 0);

   parray.append(array[8]); // 0
   TEST_ASSERT(parray.capacity() == 64);
   parray.append(array[4]); // 1
   TEST_ASSERT(parray.capacity() == 64);
   parray.append(array[3]); // 2
   parray.append(array[5]); // 3
   parray.append(array[6]); // 4
   parray.append(array[9]); // 5
   TEST_ASSERT(parray.capacity() == 64);

   TEST_ASSERT(parray.size() == 6);
   TEST_ASSERT(parray[0] == array[8]);
   TEST_ASSERT(parray[1] == array[4]);
   TEST_ASSERT(parray[2] == array[3]);
   TEST_ASSERT(parray[3] == array[5]);
   TEST_ASSERT(parray[4] == array[6]);
   TEST_ASSERT(parray[5] == array[9]);

} 

void APArrayTest::testModify()
{
   printMethod(TEST_FUNC);

   parray.reserve(2);
   parray.append(array[8]); // 0
   parray.append(array[4]); // 1
   parray.append(array[3]); // 2
   parray.append(array[5]); // 3
   parray.append(array[6]); // 4
   parray.append(array[9]); // 5
   TEST_ASSERT(parray.size() == 6);

   try {
      parray[1] = 13;
   } catch (Exception e) {
      TEST_ASSERT(0);
   }

   TEST_ASSERT(parray[1] == 13);
   TEST_ASSERT(array[4]  == 13);
   TEST_ASSERT(&parray[1]  == &array[4]);
   TEST_ASSERT(parray[2] == array[3]);
   TEST_ASSERT(&parray[2] == &array[3]);
   TEST_ASSERT(parray[3] == array[5]);
   TEST_ASSERT(&parray[3] == &array[5]);

} 

void APArrayTest::testIterator()
{
   printMethod(TEST_FUNC);

   parray.begin(iterator);
   TEST_ASSERT(iterator.isEnd());
   TEST_ASSERT(!iterator.notEnd());

   parray.reserve(2);
   parray.append(array[8]);
   parray.append(array[4]);
   parray.append(array[3]);
   parray.append(array[5]);
   TEST_ASSERT(parray.size() == 4);

   parray.begin(iterator);

   TEST_ASSERT(!iterator.isEnd());
   TEST_ASSERT(iterator.notEnd());
   TEST_ASSERT(*iterator == array[8]);
   TEST_ASSERT(iterator.get() == &array[8]);
   TEST_ASSERT(&(*iterator) == &array[8]);
   ++iterator;

   TEST_ASSERT(!iterator.isEnd());
   TEST_ASSERT(iterator.notEnd());
   TEST_ASSERT(*iterator == array[4]);
   TEST_ASSERT(iterator.get() == &array[4]);
   TEST_ASSERT(&(*iterator) == &array[4]);
   ++iterator;

   TEST_ASSERT(!iterator.isEnd());
   TEST_ASSERT(iterator.notEnd());
   TEST_ASSERT(*iterator == array[3]);
   ++iterator;

   TEST_ASSERT(!iterator.isEnd());
   TEST_ASSERT(iterator.notEnd());
   TEST_ASSERT(*iterator == array[5]);
   ++iterator;

   TEST_ASSERT(iterator.isEnd());
   TEST_ASSERT(!iterator.notEnd());

}

TEST_BEGIN(APArrayTest)
TEST_ADD(APArrayTest, testAppend)
TEST_ADD(APArrayTest, testAppendEmpty)
TEST_ADD(APArrayTest, testModify)
TEST_ADD(APArrayTest, testIterator)
TEST_END(APArrayTest)

#endif
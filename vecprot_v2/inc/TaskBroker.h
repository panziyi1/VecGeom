#ifndef GEANT_TASKBROKER
#define GEANT_TASKBROKER

class GeantBasket;

class TaskBroker
{
protected:
   struct TaskData {};
public:
   typedef TaskData *Stream;
   virtual Stream GetNextStream() = 0;

   virtual void runTask(int threadid, GeantBasket *basket) = 0;

   virtual Stream launchTask(bool wait = false) = 0;
   virtual void waitForTasks() = 0;
   
   virtual long GetTotalWork() = 0;
};


#endif // GEANT_TASKBROKER


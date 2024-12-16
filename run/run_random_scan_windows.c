#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/wait.h>

int main() {
    int num_instances = 15;  // Number of instances
    char program_path[] = "M:/Users/linusong/Documents/DMDiscovery/micromegas_6.1.15/I2HDM_DM/5-D_scan.c";

    for (int i = 0; i < num_instances; i++) {
        pid_t pid = fork();
        if (pid == 0) {
            // Child process
            char instance_id[10];
            snprintf(instance_id, sizeof(instance_id), "%d", i); // Convert instance ID to string
            char *args[] = { program_path, instance_id, NULL };  // Pass instance ID as argument

            printf("Launching instance %d with ID: %s\n", i, instance_id); // Debugging info
            execvp(args[0], args);  // Execute the program with arguments
            perror("execvp failed");  // If execvp fails
            exit(EXIT_FAILURE);
        } else if (pid > 0) {
            // Parent process: add a small delay
            usleep(50000);  // Delay of 50 milliseconds
        } else {
            perror("fork failed");
            exit(EXIT_FAILURE);
        }
    }

    // Parent process waits for all child processes to complete
    while (wait(NULL) > 0);
    return 0;
}


#include "netflow.hpp"
#include <iostream>
#define appendPref(course) \
    instructors.back().addPreference((course))
int main()
{
    /**
     * I use real instructors' names for fun, but do note that the
     * preference lists are not necessarily accurate.
     * I picked them in a way that leads to an interesting result
     * in which one of Butner or Posnett cannot be assigned
     * as many courses as he desires. (That is, if Posnett is
     * assigned ECS 154A, then Butner cannot be assigned two courses.)
     */
    std::vector<Instructor> instructors;
    instructors.emplace_back("Frid", 2);
    appendPref("ECS 120");
    appendPref("ECS 122A");
    appendPref("ECS 36C");
    instructors.emplace_back("Franklin", 1);
    appendPref("ECS 120");
    instructors.emplace_back("Porquet", 1);
    appendPref("ECS 158");
    appendPref("ECS 36C");
    instructors.emplace_back("Matloff", 1);
    appendPref("ECS 132");
    appendPref("ECS 158");
    instructors.emplace_back("Eiselt", 1);
    appendPref("ECS 32B");
    instructors.emplace_back("Butner", 1);
    appendPref("ECS 32B");
    std::vector<std::string> courses{
        "ECS 32B",
        "ECS 36C",
        "ECS 120",
        "ECS 122A",
        "ECS 132",
        "ECS 158"
    };
    assignCourses(instructors, courses);
    for (const Instructor& instructor : instructors)
    {
        std::cout << "Prof. " << instructor.lastName
                  << " is assigned to teach:";
        for (const std::string& course : instructor.assignedCourses)
            std::cout << ' ' << course;
        std::cout << std::endl;
    }
}

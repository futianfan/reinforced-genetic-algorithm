"""
An object for auto-detecting and creating jobs with the proper templates.
"""

# You'll need to import the base class first

def get_all_subclasses(base_class):
    """
    Method for getting all child classes from a parent object. Taken from:
    http://stackoverflow.com/questions/3862310/how-can-i-find-all-subclasses-of-a-class-given-its-name

    Inputs:
    :param class base_class: The parent class which we are looking towards.

    Returns
    :returns: class all_subclasses: A list of classes representing the child
        classes of the base_class
    """

    all_subclasses = []

    for subclass in base_class.__subclasses__():

        all_subclasses.append(subclass)
        all_subclasses.extend(get_all_subclasses(subclass))

    return all_subclasses

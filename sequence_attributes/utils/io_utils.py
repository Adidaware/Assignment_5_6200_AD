"""
Module to start implementing for IO
"""
import sys


class FileHandler:
    """
    This class allows us to have a more robust way of handling Exceptions, all in one place thereby conforming to DRY.
    DRY (Don't Repeat Yourself) principle involves evaluating code reuse, maintainability, and the likelihood of
    repeating code patterns for similar tasks.
    """
    def __init__(self, file, mode='r', encoding="utf-8"):
        self.file = file
        self.mode = mode
        self.file_obj = None
        self.encoding = encoding

    # This approach ensures that the file is automatically closed when the block inside the with statement is exited,
    # even if exceptions are raised, thereby making resource management more robust and error handling cleaner using the
    # context manger protocol
    # The context manager object results from evaluating the expression after with.
    # In other words, expression must return an object that implements the context management protocol.
    # This protocol consists of two special methods:

    # __enter__() is called by the with statement to enter the runtime context.
    # __exit__() is called when the execution leaves the with code block.
    def __enter__(self):
        """
        In this case, .__enter__(), typically provides the setup code.
        """
        try:
            self.file_obj = open(self.file, self.mode, encoding=self.encoding)
            return self.file_obj
        except OSError as err:
            print(f"{err}\nOSError: Could not open the file: {self.file} for mode '{self.mode}'", file=sys.stderr)
            raise err
        except ValueError as err:
            print(f"{err}\nValueError: Invalid mode '{self.mode}' for file: {self.file}", file=sys.stderr)
            raise err
        except TypeError as err:
            print(f"{err}\nTypeError: Invalid type '{self.mode}' for file: {self.file}", file=sys.stderr)
            raise err

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Once the with code block finishes, .__exit__() gets called. This method typically provides the teardown
        logic or cleanup code, such as calling .close() on an open file object. That’s why the with statement is so
        useful. It makes properly acquiring and releasing resources a breeze.
        """
        self.file_obj.close()

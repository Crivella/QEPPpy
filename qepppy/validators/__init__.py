# pylint: disable=unused-argument

def flatten_iter(iterable):
    try:
        iter(iterable)
    except TypeError:
        return [iterable,]

    if isinstance(iterable, str):
        return [iterable,]

    res = []
    for i in iterable:
        res += flatten_iter(i)
    return res

def check_shape(shape):
    def func(instance, attribute, value):
        # print(instance, attribute, value, sep='\n----------------\n')
        if value is None:
            return
        if len(value) == 0:
            return

        app = (s if isinstance(s, int) else getattr(instance, s) for s in shape)

        for i,s1 in enumerate(app):
            if s1 == -1:
                continue
            if isinstance(s1, int):
                s2 = value.shape[i]
            elif isinstance(s1, str):
                s2 = getattr(instance, s1)
            else:
                raise NotImplementedError(f"Shape must be a tuple of int or str. '{s1}' is {type(s1)}")
            if s1 != s2:
                raise ValueError(f'Shape must be {app}.')

    return func

def check_allowed(allowed):
    def func(instance, attribute, value):
        if value is None:
            return
        if len(value) == 0:
            return

        for val in flatten_iter(value):
            if val not in allowed:
                raise ValueError(f'Value must be one of {allowed}.')

    return func

def post_setter(key, func):
    def wrapper(instance, attribute, value):
        if value is None:
            return
        setattr(instance, key, func(value))

    return wrapper

def converter_none(func):
    def wrapper(value):
        # print('CONVERTING', type(value), 'to ...')
        if value is None:
            return
        # print(type(func(value)))
        return func(value)

    return wrapper

#include <Python.h>
#include <arrayobject.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "pyeos.h"

typedef struct {
    PyObject_HEAD;
    t_eos_table *eos_table;
} pyEos;

/* forward declarations */
void pyEosDealloc(PyObject* self);
static PyObject* pyEos_str(pyEos* self);
PyObject* pyEos_egiven(pyEos* self, PyObject* args);
PyObject* pyEos_tgiven(pyEos* self, PyObject* args);

static PyMethodDef pyEosMethods[] =
        {
                { "egiven", (getattrofunc)pyEos_egiven, METH_VARARGS, nullptr},
                { "tgiven", (getattrofunc)pyEos_tgiven, METH_VARARGS, nullptr},
                { nullptr, nullptr, 0, nullptr }
        };


static PyTypeObject pyEosType = {
        PyVarObject_HEAD_INIT(nullptr, 0)
        "pyeos",               /* tp_name */
        sizeof(pyEos),         /* tp_basicsize */
        0,                         /* tp_itemsize */
        pyEosDealloc,                         /* tp_dealloc */
        0,                         /* tp_print */
        nullptr,                         /* tp_getattr */
        nullptr,                         /* tp_setattr */
        nullptr,                         /* tp_reserved */
        nullptr,                         /* tp_repr */
        nullptr,                         /* tp_as_number */
        nullptr,                         /* tp_as_sequence */
        nullptr,                         /* tp_as_mapping */
        nullptr,                         /* tp_hash  */
        nullptr,                         /* tp_call */
        (reprfunc)pyEos_str,                         /* tp_str */
        nullptr,                         /* tp_getattro */
        nullptr,                         /* tp_setattro */
        nullptr,                         /* tp_as_buffer */
        Py_TPFLAGS_DEFAULT,        /* tp_flags */
        "old LEAFS EOS",            /* tp_doc */
        nullptr,                         /* tp_traverse */
        nullptr,                         /* tp_clear */
        nullptr,                         /* tp_richcompare */
        0,                         /* tp_weaklistoffset */
        nullptr,                         /* tp_iter */
        nullptr,                         /* tp_iternext */
        pyEosMethods,          /* tp_methods */
        nullptr,                         /* tp_members */
        nullptr,                         /* tp_getset */
        nullptr,                         /* tp_base */
        nullptr,                         /* tp_dict */
        nullptr,                         /* tp_descr_get */
        nullptr,                         /* tp_descr_set */
        0,                         /* tp_dictoffset */
        nullptr,                         /* tp_init */
        nullptr,                         /* tp_alloc */
        nullptr,                         /* tp_new */
};

void pyEosDealloc(PyObject* self) {
    eos_deinit(((pyEos*)self)->eos_table);
    PyObject_Del(self);
}

PyObject* pyEos_str(pyEos* self) {
    return PyUnicode_FromString(self->eos_table->datafile);
}

PyObject* pyEos_egiven(pyEos* self, PyObject* args) {
    PyArrayObject *pyXnuc;
    double rho, e, temp, p, dpdr;

    if (!PyArg_ParseTuple(args, "dO!d:eos.egiven(rho, xnuc, e)", &rho, &PyArray_Type, &pyXnuc, &e)) {
        return nullptr;
    }

    if (PyArray_NDIM(pyXnuc) != 1 || PyArray_DIMS(pyXnuc)[0] != self->eos_table->nspecies) {
        printf("%d %d\n", PyArray_NDIM(pyXnuc), self->eos_table->nspecies);
        PyErr_SetString(PyExc_ValueError, "xnuc has to be a one-dimensional array containing the abundances of all species given in the speciesfile.");
        return nullptr;
    }

    eos_calc_egiven(self->eos_table, rho, (double*)PyArray_DATA(pyXnuc), e, &temp, &p, &dpdr);

    return Py_BuildValue("ddd", temp, p, dpdr);
}

PyObject* pyEos_tgiven(pyEos* self, PyObject* args) {
    PyArrayObject *pyXnuc;
    double rho, temp, e, dedt, p;

    if (!PyArg_ParseTuple(args, "dO!d:eos.tgiven(rho, xnuc, temp)", &rho, &PyArray_Type, &pyXnuc, &temp)) {
        return nullptr;
    }

    if (PyArray_NDIM(pyXnuc) != 1 || PyArray_DIMS(pyXnuc)[0] != self->eos_table->nspecies) {
        PyErr_SetString(PyExc_ValueError, "xnuc has to be a one-dimensional array containing the abundances of all species given in the speciesfile.");
        return nullptr;
    }

    eos_calc_tgiven(self->eos_table, rho, (double*)PyArray_DATA(pyXnuc), temp, &e, &dedt, &p);

    return Py_BuildValue("ddd", e, dedt, p);
}

PyObject* PyGetEosObject(t_eos_table* eos_table) {
    pyEos* eos;

    if (eos_table) {
        eos = (pyEos*)PyObject_New(pyEos, &pyEosType);
        eos->eos_table = eos_table;
        return (PyObject*)eos;
    } else {
        Py_RETURN_NONE;
    }
}

int pyConvertEos(PyObject* object, t_eos_table** eos_table) {
    if (strcmp(object->ob_type->tp_name, "pyeos") != 0) {
        PyErr_BadArgument();
        return 0;
    }

    *eos_table = ((pyEos*)object)->eos_table;
    return 1;
}

/* pyEos module */

PyObject* loadeos(PyObject *self, PyObject *args) {
    char *datafile, *speciesfile;
    t_eos_table* eos_table;
    int oldformat, ret;

    oldformat = 0;
    if (!PyArg_ParseTuple(args, "ss|i:loadeos(datafile, speciesfile, [oldformat])", &datafile, &speciesfile, &oldformat)) {
        return nullptr;
    }

    eos_table = (t_eos_table*)malloc(sizeof(t_eos_table));
    if (eos_table == nullptr) {
        PyErr_SetFromErrno(PyExc_MemoryError);
        return nullptr;
    }
    if (!oldformat) {
        ret = eos_init(eos_table, datafile, speciesfile);
    } else {
        ret = eos_init_old(eos_table, datafile, speciesfile);
    }

    if (ret != 0) {
        free(eos_table);
        PyErr_SetString(PyExc_RuntimeError, "error in EoS initialization");
        return nullptr;
    }

    return PyGetEosObject(eos_table);
}

static PyMethodDef pyeosMethods[] = {
        { "loadeos", loadeos, METH_VARARGS, "" },
        { nullptr, nullptr, 0, nullptr }
};

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "pyeos", /* m_name */
        nullptr,      /* m_doc */
        -1,                  /* m_size */
        pyeosMethods,        /* m_methods */
        nullptr,                /* m_reload */
        nullptr,                /* m_traverse */
        nullptr,                /* m_clear */
        nullptr,                /* m_free */
};

PyMODINIT_FUNC PyInit_pyeos(void)
{
    import_array()
    pyEosType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&pyEosType) < 0)
        return nullptr;
    return PyModule_Create(&moduledef);
}
